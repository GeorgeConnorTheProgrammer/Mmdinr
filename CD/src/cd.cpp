#include <iostream>
#include <iomanip>
#include <numeric>
#include <array>

#define _USE_MATH_DEFINES
#include <cmath>

#include "cd.hpp"
#include "signedarray.hpp"

//-----------------------------------------------------------------
// Simulation Functions
//-----------------------------------------------------------------

void CDState::PrintReactionRates() {
  std::cerr << "\n[\n";
  for (int i = -MAX_SIZE; i < MAX_SIZE; ++i) {
    for (int j = -MAX_SIZE; j < MAX_SIZE; ++j) {
      std::cerr << std::setfill(' ') << std::setw(10) << reaction_rates[i][j] << ", ";
    }
    std::cerr << "\n";
  }
  std::cerr << "]\n";
}

void CDState::Init()
{
  // Based on the example case in section 4.4 of the Kohnert paper

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

    double E_b = 1.73 - 2.59 * (std::pow(i, 2/3) - std::pow(i-1, 2/3)); // TODO - this is only really correct for vacancies
    int dissociation_direction = i > 0 ? -1 : 1;
    dissociation_rates[i][i + dissociation_direction] = ((reaction_rates[i][i + dissociation_direction] + reaction_rates[i + dissociation_direction][i]) / atomic_volume) * std::exp(-E_b / (k * T));
  }

  species[1].g = 1000;
  species[1].r_s = std::pow(10, 3);
  species[1].K = 4 * M_PI * species[1].r_s * species[1].D;

  species[-1].g = 0.01;
  species[-1].r_s = std::pow(10, 3);
  species[-1].K = 4 * M_PI * species[1].r_s * species[1].D;
  
  prev_species.set(species);
}

double CDState::GetReactionRate(int i) {
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
    double dissociation_rate = 0.0; // TODO

    return combination_rate + dissociation_rate;
}

void CDState::Step(double dt)
{
  for (int i = -MAX_SIZE; i <= MAX_SIZE; ++i) // Compute the change in concentration for each cluster species
  {
    if (i == 0) continue;

    Species& s = species[i];

    double R = GetReactionRate(i);
    double sink_loss = s.K * C_s * s.C;
    double dC = dt * (s.g + R - sink_loss);

    s.C += dC;
  }

  prev_species.set(species);
}