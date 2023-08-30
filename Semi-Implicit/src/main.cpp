#include <fstream>

#include "Navier2.hpp"
#include "Preconditioners.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  const double T = 2.0;

  double deltat = 0.05;

  Navier problem(degree_velocity, degree_pressure, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}