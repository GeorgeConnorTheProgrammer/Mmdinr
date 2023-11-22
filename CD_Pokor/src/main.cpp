#include <iostream>

#include "cd.hpp"

void runCD(double dt, double total_time)
{
  CDState cd;

  cd.Init();

  for (double t = 0; t < total_time; t += dt)
  {
    cd.Step(dt);
  }
}

int main(int argc, char** argv)
{
  if (argc < 3) 
  {
    std::cout << "Too few args. Usage: test [dt] [steps]" << std::endl;
    return 1;
  }

  double dt = atof(argv[1]);
  double total_time = atof(argv[2]);

  runCD(dt, total_time);
  return 0;
}