#include <iostream>
#include <numeric>
#include <array>

namespace CD
{
}

// Physical parameters

double C_i = 0.0; // Concentration of interstitials
double C_v = 0.0; // Concentration of vacancies

double g_i = 0.1; // Generation of interstitials due to radiation
double g_v = 0.1; // Generation of vacancies due to radiation

double D_i = 0.4; // Interstitial diffusion coefficient
double D_v = 0.6; // Vacancy diffusion coefficient

double k_iv = 0.5; // Recombination rate constant

// Interstitial sink partial strengths
std::array<double, 3> i_sinks {
  0.01, 0.12, 0.28
};

// Vacancy sink partial strengths
std::array<double, 3> v_sinks {
  0.01, 0.8, 0.18
};

double k_i_2 = 0.0; // Total sink strength for removal of interstitials
double k_v_2 = 0.0; // Total sink strength for removal of vacancies

// Additional simulation state

double t = 0.0;

// Run the simulation

void runModel(double dt, int steps) 
{
  k_i_2 = std::accumulate(i_sinks.begin(), i_sinks.end(), 0.0, std::plus<double>());
  k_v_2 = std::accumulate(v_sinks.begin(), v_sinks.end(), 0.0, std::plus<double>());

  std::cout << "Ci, Cv" << std::endl;

  for (int i = 0; i < steps; ++i) 
  {
    double dC_i = (g_i - k_iv * C_i * C_v - D_i * k_i_2 * C_i) * dt;
    double dC_v = (g_v - k_iv * C_i * C_v - D_v * k_v_2 * C_v) * dt;

    C_i += dC_i;
    C_v += dC_v;

    std::cout << "" << C_i << "," << C_v << std::endl;
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
  int steps = atoi(argv[2]);

  runModel(dt, steps);
  return 0;
}