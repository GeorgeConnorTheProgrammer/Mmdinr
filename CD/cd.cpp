#include <iostream>
#include <iomanip>
#include <numeric>
#include <array>

#define _USE_MATH_DEFINES
#include <cmath>

#include "signedarray.hpp"


//-----------------------------------------------------------------
// Simulation Classes
//-----------------------------------------------------------------

struct Species
{
  double C = 0.0; // Concentration of this species
  double g = 0.0; // Generation of this cluster species in collision cascades
  double D = 0.0; // Diffusion constant
  double K = 0.0; // Sink strength
  
  double r = 0.0; // Reaction radius
  double r_s = 0.0; // Reaction radius with sinks
};




//-----------------------------------------------------------------
// Simulation State
//-----------------------------------------------------------------

constexpr int MAX_SIZE = 20;
SignedArray<Species, MAX_SIZE> species{};
SignedArray<Species, MAX_SIZE> prev_species{};

SignedArray<SignedArray<double, MAX_SIZE>, MAX_SIZE> reaction_rates;
SignedArray<SignedArray<double, MAX_SIZE>, MAX_SIZE> dissociation_rates;


//-----------------------------------------------------------------
// Simulation Functions
//-----------------------------------------------------------------

void print_reaction_rates() {
  std::cerr << "\n[\n";
  for (int i = -MAX_SIZE; i < MAX_SIZE; ++i) {
    for (int j = -MAX_SIZE; j < MAX_SIZE; ++j) {
      std::cerr << std::setfill(' ') << std::setw(10) << reaction_rates[i][j] << ", ";
    }
    std::cerr << "\n";
  }
  std::cerr << "]\n";
}

void initModel()
{
  // Based on the example case in section 4.4 of the Kohnert paper

  for (int i = -reaction_rates.size(); i < reaction_rates.size(); ++i) {
    for (int j = -reaction_rates.size(); j < reaction_rates.size(); ++j) {
      reaction_rates[i][j] = -42.0;
    }
  }

  const double T = 300; // Temperature in Kelvin
  const double atomic_volume = 0.0118; // Atomic volume in nm^3

  const double E_mv = 0.67; //Migration energy of point vacancies in eV
  const double E_mi = 0.34; //Migration energy of point interstitials in eV
  const double k = 8.6173 * std::pow(10,-5); //eV K^-1 k is the Boltzmann constant

  species[1].D = std::pow(10, 11) * std::exp(-E_mi / (k * T));
  species[-1].D = std::pow(10, 11) * std::exp(-E_mv / (k * T));

  for (int i = -MAX_SIZE; i < MAX_SIZE; ++i) {
    if (i == 0) continue;

    species[i].r = std::cbrt(3 * std::abs(i) * atomic_volume / (4 * M_PI)); // TODO - Should we really assume all clusters are spherical?

    for (int j = -MAX_SIZE; j < MAX_SIZE; ++j) { // Calculate reaction rate coefficients
      if (j == 0) continue;

      if (i + j < MAX_SIZE && i + j > -MAX_SIZE) {
        const double r_ij = species[i].r + species[j].r;
        reaction_rates[i][j] = 4 * M_PI * r_ij * species[i].D;
        reaction_rates[j][i] = 4 * M_PI * r_ij * species[j].D;
      }
    }

    dissociation_rates[i][i-1] = 0.0; // TODO
    dissociation_rates[i][i-1] = 0.0;
  }

  species[1].g = 0.001;
  species[1].r_s = std::pow(10, 3);
  species[1].K = 4 * M_PI * species[1].r_s * species[1].D;

  species[-1].g = 0.001;
  species[-1].r_s = std::pow(10, 3);
  species[-1].K = 4 * M_PI * species[1].r_s * species[1].D;
  
  prev_species.set(species);
}

double calculate_reaction_rate(int i) {
    double j_plus_k_equals_i = 0.0; // Run through all reactions that can create i
    for (int j = -MAX_SIZE; j <= MAX_SIZE; ++j)
    {
      if (j == 0 || j == i) continue;
      double k = i - j;
      if (k < -MAX_SIZE || k > MAX_SIZE) continue;
      j_plus_k_equals_i += reaction_rates[j][k] * species[j].C * species[k].C;
    }

    double i_plus_j_equals_k = 0.0; // Run through all reactions that can be created using i
    for (int j = -MAX_SIZE; j <= MAX_SIZE; ++j)
    {
      if (j == 0) continue;
      double k = i + j;
      if (k < -MAX_SIZE || k > MAX_SIZE) continue;
      i_plus_j_equals_k -= reaction_rates[i][j] * species[i].C * species[j].C;
    }
    double combination_rate = j_plus_k_equals_i - i_plus_j_equals_k;
    
    double dissociation_rate = 0.0;

    return combination_rate + dissociation_rate;
}

void runStep(double dt)
{
  for (int i = -MAX_SIZE; i <= MAX_SIZE; ++i) // Compute the change in concentration for each cluster species
  {
    if (i == 0) continue;

    Species& s = species[i];

    double R = calculate_reaction_rate(i);
    double dC = dt * (s.g + R);

    s.C += dC;
  }

  prev_species.set(species);
}

void runModel(double dt, double total_time)
{
  initModel();
  print_reaction_rates();

  std::cout << "t";
  for (int i = -MAX_SIZE; i <= MAX_SIZE; ++i)
  {
    if (i == 0) continue;
    std::cout << ", C_" << i; 
  }
  std::cout << "\n";

  for (double t = 0; t < total_time; t += dt)
  {
    runStep(dt);
  }

  std::cout << total_time;
  for (int i = -MAX_SIZE; i <= MAX_SIZE; ++i) {
    if (i == 0) continue;
    std::cout << ", " << std::log(species[i].C + 1);
  }
  std::cout << "\n";
}




//-----------------------------------------------------------------
// Main
//-----------------------------------------------------------------

int main(int argc, char** argv)
{
  if (argc < 3) 
  {
    std::cout << "Too few args. Usage: test [dt] [steps]" << std::endl;
    return 1;
  }

  double dt = atof(argv[1]);
  double total_time = atof(argv[2]);

  runModel(dt, total_time);
  return 0;
}
