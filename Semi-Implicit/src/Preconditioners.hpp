#ifndef PRECONDITIONERS_HPP
#define PRECONDITIONERS_HPP

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

using namespace dealii;

/// The SIMPLE preconditioner.
// Psimple = P1 * P2 = [F 0; B -S] * [I D^-1*B_T; 0 alpha*I] with S = B * D^-1 * B_T the approx. Schur complement
class PreconditionSIMPLE
{
public:
    // Initialize the preconditioner, given:
    // - the F matrix = M/deltat + A + C(u_n)
    // - the B matrix
    // - the B_T matrix
    // For simplicity, just pass the system matrix and extract the blocks.
    // - the alpha parameter (relaxation parameter to dump the pressure update)
    void
    initialize(const double &alpha_,
               const TrilinosWrappers::BlockSparseMatrix &system_matrix)
    {
        alpha = &alpha_;
        F = &system_matrix.block(0, 0);
        B = &system_matrix.block(1, 0);
        B_T = &system_matrix.block(0, 1);

        preconditioner_F.initialize(system_matrix.block(0, 0));
        preconditioner_B.initialize(system_matrix.block(1, 0));
        preconditioner_B_T.initialize(system_matrix.block(0, 1));
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
        // Effect of P1
        {

            // dst_0 = F^-1 * src_0
            SolverControl solver_control_F(1000,
                                           1e-2 * src.block(0).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F(
                solver_control_F);
            solver_gmres_F.solve(*F,
                                 dst.block(0),
                                 src.block(0),
                                 preconditioner_F);

            // tmp = -B * dst_0 + src_1
            tmp.reinit(src.block(1));
            B->vmult(tmp, dst.block(0));
            tmp.sadd(-1.0, src.block(1));

            // dst_1 = - B_T^-1 * D * B^-1 * tmp

            // 1. Solve B * dst_1 = tmp
            SolverControl solver_control_S(1000,
                                           1e-2 * src.block(1).l2_norm());
            SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_S(
                solver_control_S);
            solver_gmres_S.solve(*B,
                                 dst.block(1),
                                 tmp,
                                 preconditioner_B);

            // 2. extract D and tmp = D * dst_1
            D.reinit(src.block(0));
            for (unsigned int i = 0; i < D.size(); ++i)
                D(i) = F->diag_element(i);
            tmp = D * dst.block(1);
            // D->vmult(tmp, dst.block(1));

            // 3. Solve B_T * dst_1 = tmp
            solver_gmres_S.solve(*B_T,
                                 dst.block(1),
                                 tmp,
                                 preconditioner_B_T);

            // Finally, dst_1 = - dst_1
            dst.block(1) *= -1.0;
        }

        // Effect of P2
        {

            // dst_1 = 1/alpha * dst_1
            dst.block(1) *= 1.0 / *alpha;

            // tmp = B_T * dst_1
            tmp.reinit(src.block(0));
            B_T->vmult(tmp, dst.block(1));

            // tmp = D^-1 * tmp

            // Invert the diagonal
            for (unsigned int i = 0; i < D.size(); ++i)
            {
                if (D(i) != 0.0)
                    D(i) = 1.0 / D(i);
            }

            tmp = D * tmp;

            // dst_0 = dst_0 - tmp
            dst.block(0) -= tmp;
        }
    }

protected:
    // Alpha parameter.
    const double *alpha;

    // F matrix.
    const TrilinosWrappers::SparseMatrix *F;

    // the diagonal part of F
    mutable TrilinosWrappers::MPI::Vector D;

    // Preconditioner used for the F block.
    TrilinosWrappers::PreconditionILU preconditioner_F;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Preconditioner used for B.
    TrilinosWrappers::PreconditionILU preconditioner_B;

    // B_T matrix.
    const TrilinosWrappers::SparseMatrix *B_T;

    // Preconditioner used for B_T.
    TrilinosWrappers::PreconditionILU preconditioner_B_T;

    // Temporary vectors.
    mutable TrilinosWrappers::MPI::Vector tmp;
};

#endif