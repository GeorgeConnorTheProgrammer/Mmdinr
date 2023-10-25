#include <iostream>
#include <numeric>
#include <array>

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
std::array<double, NUM_SPECIES * NUM_SPECIES> dissociation_rates{};

void populateSpecies()
{
  for (int i = 0; i < NUM_SPECIES; ++i) {
    species[i].C = 0.0;
    species[i].div_flux = 0.1;
    species[i].g = 0.2;
    species[i].D = 0.5;
    species[i].k_2 = 0.01;
  }

  // TODO - very important, give these values that are at least mass-conservative
  for (int i = 0; i < NUM_SPECIES * NUM_SPECIES; ++i) {
    reaction_rates[i] = 0.1;
    dissociation_rates[i] = 0.0;
  }
}

void runStep(double dt)
{
  std::copy(species.begin(), species.end(), out_species.begin());

  for (int i = 0; i < species.size(); ++i) // Compute the change in concentration for each cluster species
  {
    Species& s = species[i];

    double combination_rate = 0.0; // Rate of change of concentration of species i due to combination of other clusters
    double dissociation_rate = 0.0; // Rate of change of concentration of species i due to dissociation

    for (int i = 0; i < species.size(); ++i) {
      for (int j = 0; j < species.size(); ++j) {
        if (i == j) continue;

        combination_rate += reaction_rates[NUM_SPECIES * j + i] * species[i].C * species[j].C; // Species j becomes species i
        combination_rate -= reaction_rates[NUM_SPECIES * i + j] * species[j].C * species[i].C; // Species i becomes species j
        
	dissociation_rate += dissociation_rates[NUM_SPECIES * j + i] * species[j]; // Species j becomes species i
        dissociation_rate -= dissociation_rates[NUM_SPECIES * i + j] * species[i]; // Species i becomes species j
      }
    }


    double R = combination_rate + dissociation_rate;
    double dC = dt * (-s.div_flux + s.g + R - s.D * s.k_2 * s.C); // TODO - investigate ways to deal with "stiff equations"

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
