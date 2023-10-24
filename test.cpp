#include <iostream>
#include <numeric>
#include <array>
#include <cmath>


namespace CD
{
  struct Species
  {
    double C = 0.0; // Concentration of this species
    double div_flux = 0.0; // Divergence of flux. Deserves scrutiny.
    double g = 0.0; // Generation of this cluster species in collision cascades
    double D = 0.0; // Diffusion constant
    double k_2 = 0.0; // Sink strength
  };

  constexpr size_t NUM_SPECIES = 2;
  std::array<Species, NUM_SPECIES> species{};
  std::array<Species, NUM_SPECIES> out_species{};

  std::array<double, NUM_SPECIES * NUM_SPECIES> reaction_rates{};

  void populateSpecies()
  {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      species[i].C = 0.0;
      species[i].div_flux = 0.1;
      species[i].g = 0.2;
      species[i].D = 0.5;
      species[i].k_2 = 0.01;
    }

    for (int i = 0; i < reaction_rates.size(); ++i) {
      reaction_rates[i] = 0.1;
    }
  }

  void runStep(double dt)
  {
    std::copy(species.begin(), species.end(), out_species.begin());

    for (int i = 0; i < species.size(); ++i) 
    {
      Species& s = species[i];

      // TODO
      double R_gain_by_combination = 0.0;
      double R_loss_by_combination = 0.0;
      double R_gain_by_dissociation = 0.0;
      double R_loss_by_dissociation = 0.0;

      for (int i = 0; i < species.size(); ++i) {
        for (int j = 0; j < species.size(); ++j) {
          if (i == j) continue;

          R_gain_by_combination += dt * reaction_rates[NUM_SPECIES * j + i] * species[i].C + species[j].C;
          
          R_loss_by_combination += 0.0;
          R_gain_by_dissociation += 0.0;
          R_loss_by_dissociation += 0.0;
        }
      }


      double R = R_gain_by_combination + R_loss_by_combination
               + R_gain_by_dissociation + R_loss_by_dissociation;

      double dC = dt * (-s.div_flux + s.g + R - s.D * s.k_2 * s.C);

      out_species[i].C += dC;
    }
  }

  void runModel(double dt, int steps)
  {

    populateSpecies();

    for (int i = 0; i < steps; ++i)
    {
      runStep(dt);
    }

  }
}

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
