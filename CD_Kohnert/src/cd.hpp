#pragma once

#include <cmath>

#include "signedarray.hpp"

struct Species
{
  double C = 0.0; // Concentration of this species
  double g = 0.0; // Generation of this cluster species in collision cascades
  double D = 0.0; // Diffusion constant
  double K = 0.0; // Sink strength
  
  double r = 0.0; // Reaction radius
  double r_s = 0.0; // Reaction radius with sinks
};

class CDState 
{
public:
  static constexpr double T = 300.0; // Temperature in Kelvin
  static constexpr double atomic_volume = 0.0118; // Atomic volume in nm^3
  static constexpr double C_s = 8.0 * 0.000001;

  static constexpr double E_mv = 0.67; //Migration energy of point vacancies in eV
  static constexpr double E_mi = 0.34; //Migration energy of point interstitials in eV
  static constexpr double k = 8.6173 * 0.00005; //eV K^-1 k is the Boltzmann constant

  static constexpr int MAX_SIZE = 40;
  SignedArray<Species, MAX_SIZE> species{};
  SignedArray<Species, MAX_SIZE> prev_species{};

  SignedArray<SignedArray<double, MAX_SIZE>, MAX_SIZE> reaction_rates;
  SignedArray<SignedArray<double, MAX_SIZE>, MAX_SIZE> dissociation_rates;

  void Init(); // TODO - Get input parameters through here
  void Step(double dt);

  double GetReactionRate(int i);
  void PrintReactionRates();
};