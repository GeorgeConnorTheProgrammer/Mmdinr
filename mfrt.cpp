#include <iostream>
#include <numeric>
#include <array>
#include <cmath>

namespace MFRT
{
  // Physical parameters

  double C_i = 0.0; // Concentration of interstitials
  double C_v = 0.0; // Concentration of vacancies

  double g_i = 0.1; // Generation of interstitials due to radiation
  double g_v = 0.1; // Generation of vacancies due to radiation
  
  // using Material parameters from SA304
  double D_0i = 0.010; //cm2/s
  double D_0v = 0.6; //cm2/s
  double E_mv = 1.35; //eV
  double E_mi = 0.45; //eV
  double k = 8.6173 * std::pow(10,-5); //eV K^-1 k is the Boltzmann constant
  double Temp = 300; //Temperature in Kelvin
  
  double D_i = D_0i * std::exp(-(E_mi/(k*Temp))); // Interstitial diffusion coefficient
  double D_v = D_0v * std::exp(-(E_mv/(k*Temp))); // Vacancy diffusion coefficient
  
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

  MFRT::runModel(dt, steps);

  return 0;

}
