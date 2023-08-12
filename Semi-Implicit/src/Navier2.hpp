#ifndef NAVIER2_HPP
#define NAVIER2_HPP

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
#include "Preconditioners.hpp"

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
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = alpha * get_time();
      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int component = 0) const override
    {
      if (component == 0)
        return alpha * get_time();
      else
        return 0.0;
    }

  protected:
    const double alpha = 1.0;
  };

  // The initial solution to set up the system.
  class InitialSolution : public Function<dim>
  {
  public:
    InitialSolution()
        : Function<dim>(dim + 1)
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // A class to possibly specify the boundary conditions on the Neumann boundary surfaces.
  class BoundaryNeumann : public Function<dim>
  {
  public:
    BoundaryNeumann()
        : Function<dim>(dim + 1)
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  Navier(const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_,
         const double &T_,
         const double &deltat_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), T(T_), degree_velocity(degree_velocity_), degree_pressure(degree_pressure_), deltat(deltat_), mesh(MPI_COMM_WORLD)
  {
  }

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

  // Compute the error.
  // double
  //   compute_error(const VectorTools::NormType& norm_type);

protected:
  // Assemble the system for the time step, using the semi-implicit treatment of the non linear convective term.
  void
  assemble(const double &time);

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output.
  void
  output(const unsigned int &time_step, const double &time) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  double nu = 0.1;

  // Gamma parameter
  const double alpha = 0.5;

  // Outlet pressure [Pa].
  const double p_out = 10;

  // NEW PARAMETERS FOR REYNOLDS NUMBER APPROXIMATION //////////////////////////

  // Characteristic length [m].
  const double L = 1.0; // In this case, the leght of the cylinder

  // Characteristic velocity [m/s].
  const double U = 10.0; // In this case, the inlet velocity

  // Density [kg/m3].
  const double rho = 1.225; // In this case, the density of the liquid flowing around the cylinder

  /////////////////////////////////////////////////////////////////////////////

  // Forcing term.
  ForcingTerm forcing_term;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Boundary conditions.
  // BoundaryNeumann boundary_neumann;

  // Initial solution.
  InitialSolution initial_solution;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Time step.
  const double deltat;

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
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Diagonal Inverse of F matrix, needed for preconditioning. We use a blockvector for comodity.
  TrilinosWrappers::MPI::BlockVector D_inv;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::BlockVector solution_old;
};

#endif