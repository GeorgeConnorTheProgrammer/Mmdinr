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