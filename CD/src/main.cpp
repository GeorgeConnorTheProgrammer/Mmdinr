#include <iostream>

#include "cd.hpp"

void runCD(double dt, double total_time)
{
  CDState cd;

  cd.Init();
  //cd.PrintReactionRates();

  std::cout << "t";
  for (int i = -cd.MAX_SIZE; i <= cd.MAX_SIZE; ++i)
  {
    if (i == 0) continue;
    std::cout << ", C_" << i; 
  }
  std::cout << "\n";

  for (double t = 0; t < total_time; t += dt)
  {
    std::cout << t;
    for (int i = -cd.MAX_SIZE; i <= cd.MAX_SIZE; ++i) {
      if (i == 0) continue;
      std::cout << ", " << std::log(cd.species[i].C + 1);
    }
    std::cout << "\n";

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