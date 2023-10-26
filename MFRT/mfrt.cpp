#include <iostream>
#include <fstream>
#include <numeric>
#include <array>
#include <cmath>

#include "../vendor/nlohmann/json.hpp"

namespace MFRT
{
  // Physical parameters

  double C_i = 0.0; // Concentration of interstitials
  double C_v = 0.0; // Concentration of vacancies
 
  // using Material parameters from SA304
  double D_0i = 0.001; //cm2/s
  double D_0v = 0.6; //cm2/s
  double E_mv = 1.35; //eV
  double E_mi = 0.45; //eV
  double k = 8.6173 * std::pow(10,-5); //eV K^-1 k is the Boltzmann constant
	double r_iv = 7 * std::pow(10,-8); //should be 7nm but the simulation is in cm

	// We can assume r_vs = r_is = 10^-4 cm according to 10-19-23 slides.
	double r_vs = std::pow(10,-4);
	double r_is = r_vs;

  // Run the simulation
  void run_model(double dt, int steps, double Temp, int K_0_exp, int C_s_exp) 
  {
		// running variable init
		// Calculating the D_i and D_v.
		double D_i = D_0i * std::exp(-(E_mi/(k*Temp))); // Interstitial diffusion coefficient
		double D_v = D_0v * std::exp(-(E_mv/(k*Temp))); // Vacancy diffusion coefficient		
		
		double K_0 = std::pow(10,K_0_exp); //defect production rate
		double C_s = std::pow(10,C_s_exp); //Mentioned in 10-19-23 rate theory slides, needs clarification..

		double K_iv = 4.0 * M_PI * r_iv * (D_i + D_v); // vancancy-interstitial recombination rate coeff
		double K_is = 4.0 * M_PI * r_is * D_i; // interstitial-sink reaction rate coeff.
		double K_vs = 4.0 * M_PI * r_vs * D_v; // vancancy-sink reaction rate coeff.
    
		std::cout << "Ci, Cv" << std::endl;

    for (int i = 0; i < steps; ++i) 
    {
      double dC_i = (K_0 - K_iv * C_i * C_v - K_vs * C_i * C_s) * dt;
      double dC_v = (K_0 - K_iv * C_i * C_v - K_is * C_i * C_s) * dt;
			if (std::isinf(dC_i) || std::isinf(dC_v) || dC_i != dC_i || dC_v != dC_v)
			{
				std::cout << "limit reached stopping model.." << std::endl;
				break; 
			}
      C_i += dC_i;
      C_v += dC_v;

      std::cout << "" << C_i << "," << C_v << std::endl;
    }
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) 
  {
    std::cerr << "Too few args. Usage: mfrt [config file name]" << std::endl;
		return 1;
  }

  std::ifstream config_file(argv[1]);
  if (!config_file.good()) {
    std::cerr << "Could not open " << argv[1] << std::endl;
  }
  nlohmann::json config = nlohmann::json::parse(config_file);

  double time = config["total_time_seconds"];
  double dt = config["dt_seconds"];
	double Temp = config["temperature_kelvin"];
	int K_0_exp = config["K_0_exp"];
	int C_s_exp = config["C_s_exp"];

  MFRT::run_model(dt, time, Temp, K_0_exp, C_s_exp);

  return 0;
}
