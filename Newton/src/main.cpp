#include "Navier.hpp"

// Main function.
int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int degree_velocity = 2;
    const unsigned int degree_pressure = 1;

    Navier problem(degree_velocity, degree_pressure);

    problem.setup();
    problem.run();
    problem.output();

    return 0;
}