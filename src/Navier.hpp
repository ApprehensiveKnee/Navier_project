#ifndef NAVIER_HPP
#define NAVIER_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Navier-Stokes problem.
class Navier
{
public:
    // Physical dimension (1D, 2D, 3D)
    static constexpr unsigned int dim = 3;

    // Function for the forcing term.
    class ForcingTerm : public Function<dim>
    {
    public:
        virtual void
        vector_value(const Point<dim> & /*p*/,
                     Vector<double> &values) const override
        {
            for (unsigned int i = 0; i < dim - 1; ++i)
                values[i] = 0.0;

            values[dim - 1] = -g;
        }

        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int component = 0) const override
        {
            if (component == dim - 1)
                return -g;
            else
                return 0.0;
        }

    protected:
        const double g = 0.0;
    };

    // Function for inlet velocity. This actually returns an object with four
    // components (one for each velocity component, and one for the pressure), but
    // then only the first three are really used (see the component mask when
    // applying boundary conditions at the end of assembly). If we only return
    // three components, however, we may get an error message due to this function
    // being incompatible with the finite element space.
    class InletVelocity : public Function<dim>
    {
    public:
        InletVelocity()
            : Function<dim>(dim + 1)
        {
        }

        virtual void
        vector_value(const Point<dim> &p, Vector<double> &values) const override
        {
            // values[0] = -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
            values[0] = alpha;
            for (unsigned int i = 1; i < dim + 1; ++i)
                values[i] = 0.0;
        }

        virtual double
        value(const Point<dim> &p, const unsigned int component = 0) const override
        {
            if (component == 0)
                // return -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
                return alpha;
            else
                return 0.0;
        }

    protected:
        const double alpha = 10.0;
    };

    // Since we're working with block matrices, we need to make our own
    // preconditioner class. A preconditioner class can be any class that exposes
    // a vmult method that applies the inverse of the preconditioner.

    // Identity preconditioner.
    class PreconditionIdentity
    {
    public:
        // Application of the preconditioner: we just copy the input vector (src)
        // into the output vector (dst).
        void
        vmult(TrilinosWrappers::MPI::BlockVector &dst,
              const TrilinosWrappers::MPI::BlockVector &src) const
        {
            dst = src;
        }

    protected:
    };

    // Block-triangular preconditioner. (Block Shur Preconditioner, as described in https://www.dealii.org/current/doxygen/deal.II/step_57.html)
    class PreconditionBlockTriangular
    {
    public:
        // Initialize the preconditioner, given:
        // - the velocity stiffness matrix
        // - the pressure mass matrix
        // - the B matrix
        // - the viscosity
        // - the gamma parameter
        void
        initialize(const double &gamma_,
                   const double &viscosity_,
                   const TrilinosWrappers::SparseMatrix &velocity_stiffness_modified_,
                   const TrilinosWrappers::SparseMatrix &pressure_mass_,
                   const TrilinosWrappers::SparseMatrix &B_T_)
        {
            gamma = &gamma_;
            viscosity = &viscosity_;
            velocity_stiffness = &velocity_stiffness_modified_;
            pressure_mass = &pressure_mass_;
            B_T = &B_T_;

            preconditioner_velocity.initialize(velocity_stiffness_modified_);
            preconditioner_pressure.initialize(pressure_mass_);
        }

        // Application of the preconditioner.
        void
        vmult(TrilinosWrappers::MPI::BlockVector &dst,
              const TrilinosWrappers::MPI::BlockVector &src) const
        {

            // Effect lower part of the preconditioner (Shur complement)
            // Computing dst_1 = -(gamma + viscosity)*Mp^-1 * src_1
            {
                SolverControl solver_control_pressure(1000,
                                                      1e-2 * src.block(1).l2_norm());
                SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
                    solver_control_pressure);
                solver_cg_pressure.solve(*pressure_mass,
                                         dst.block(1),
                                         src.block(1),
                                         preconditioner_pressure);
                dst.block(1) *= -(*gamma + *viscosity);
            }

            // Effect of the upper part of the block-triangular preconditioner
            // Computing dest_1 = A^-1 * (src_0 - B^T * dst_1)

            {
                tmp.reinit(src.block(0));
                B_T->vmult(tmp, dst.block(1));
                tmp.sadd(-1.0, src.block(0));

                // Change the iterative solver for the velocity sub-block (A is not simmetric anymore)
                SolverControl solver_control_velocity(1000,
                                                      1e-2 * src.block(0).l2_norm());
                SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
                    solver_control_velocity);
                solver_gmres_velocity.solve(*velocity_stiffness,
                                            dst.block(0),
                                            tmp,
                                            preconditioner_velocity);
            }
        }

    protected:
        // Gamma parameter.
        const double *gamma;

        // Viscosity.
        const double *viscosity;

        // Velocity stiffness matrix.
        const TrilinosWrappers::SparseMatrix *velocity_stiffness;

        // Preconditioner used for the velocity block.
        TrilinosWrappers::PreconditionILU preconditioner_velocity;

        // Pressure mass matrix.
        const TrilinosWrappers::SparseMatrix *pressure_mass;

        // Preconditioner used for the pressure block.
        TrilinosWrappers::PreconditionILU preconditioner_pressure;

        // B matrix.
        const TrilinosWrappers::SparseMatrix *B_T;

        // Temporary vector.
        mutable TrilinosWrappers::MPI::Vector tmp;
    };

    // Constructor.
    Navier(const unsigned int &N_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_)
        : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), N(N_), degree_velocity(degree_velocity_), degree_pressure(degree_pressure_), mesh(MPI_COMM_WORLD)
    {
    }

    // Setup system.
    void
    setup();

    // Assemble system for each iteration of the Newton method. We also assemble the pressure mass matrix (needed for the
    // preconditioner).
    void
    assemble_system(const bool initial_step);

    // Solve Newton step.
    void
    solve_newton_step();

    // Solve Newton method.
    void
    solve_newton();

    // Output results.
    void
    output();

protected:
    // MPI parallel. /////////////////////////////////////////////////////////////

    // Number of MPI processes.
    const unsigned int mpi_size;

    // This MPI process.
    const unsigned int mpi_rank;

    // Parallel output stream.
    ConditionalOStream pcout;

    // Problem definition. ///////////////////////////////////////////////////////

    // Kinematic viscosity [m2/s].
    const double nu = 1;

    // Gamma parameter
    const double gamma = 1.0;

    // Outlet pressure [Pa].
    const double p_out = 10;

    // Forcing term.
    ForcingTerm forcing_term;

    // Inlet velocity.
    InletVelocity inlet_velocity;

    // Discretization. ///////////////////////////////////////////////////////////

    // Mesh refinement.
    const unsigned int N;

    // Polynomial degree used for velocity.
    const unsigned int degree_velocity;

    // Polynomial degree used for pressure.
    const unsigned int degree_pressure;

    // Mesh.
    parallel::fullydistributed::Triangulation<dim> mesh;

    // Finite element space.
    std::unique_ptr<FiniteElement<dim>> fe;

    // Quadrature formula.
    std::unique_ptr<Quadrature<dim>> quadrature;

    // Quadrature formula for face integrals.
    std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;

    // DoFs owned by current process in the velocity and pressure blocks.
    std::vector<IndexSet> block_owned_dofs;

    // DoFs relevant to the current process (including ghost DoFs).
    IndexSet locally_relevant_dofs;

    // DoFs relevant to current process in the velocity and pressure blocks.
    std::vector<IndexSet> block_relevant_dofs;

    // System matrix.
    TrilinosWrappers::BlockSparseMatrix jacobian_matrix;

    // Pressure mass matrix, needed for preconditioning. We use a block matrix for
    // convenience, but in practice we only look at the pressure-pressure block.
    TrilinosWrappers::BlockSparseMatrix pressure_mass;

    // Right-hand side vector in the linear system.
    TrilinosWrappers::MPI::BlockVector residual_vector;

    // Solutin increment (without ghost elements).
    TrilinosWrappers::MPI::BlockVector delta_owned;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::BlockVector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::BlockVector solution;
};

#endif