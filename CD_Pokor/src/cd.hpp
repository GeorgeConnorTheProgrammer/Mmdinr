#pragma once

#include <array>
#include <cmath>

class CDState 
{
public:
  static constexpr size_t NUM_CLUSTER_SIZES = 20;
  
  std::array<double, NUM_CLUSTER_SIZES> i_concentrations {0.0};
  std::array<double, NUM_CLUSTER_SIZES> v_concentrations {0.0};

  void Init(); // TODO - Get input parameters through here
  void Step(double dt);

  double C_i(int n);
  double C_v(int n);

private:
  static constexpr double k = 8.6173 * 0.00005; //eV K^-1 k is the Boltzmann constant

  static constexpr double T = 300.0; // Temperature in Kelvin
  static constexpr double r_iv = 0.00000007; // i/v reaction radius in cm
  static constexpr double D_0i = 0.001; // cm^2/s
  static constexpr double D_0v = 0.6; // cm^2/s
  static constexpr double E_mi = 0.45; // Interstitial migration energy in eV
  static constexpr double E_mv = 1.35; // Interstitial migration energy in eV

  std::array<double, NUM_CLUSTER_SIZES> i_concentrations_swap {0.0};
  std::array<double, NUM_CLUSTER_SIZES> v_concentrations_swap {0.0};

  double rho = 0.0; // Dislocation network density

  double dCi1(double dt);
  double dCv1(double dt);

  double dCi(int n, double dt);
  double dCv(int n, double dt);
  double dRho(double dt);

  double G_i(int n); // Interstitial cluster generation term
  double G_v(int n); // Vacancy cluster generation term

  double a_i(int n); 
  double b_i(int n);
  double c_i(int n);
  double a_v(int n); 
  double b_v(int n);
  double c_v(int n);

  double alpha_ii(int n);
  double alpha_vv(int n);

  double beta_iv(int n);
  double beta_vi(int n);
  double beta_ii(int n);
  double beta_vv(int n);

  double ta_gbi();
  double ta_gbv();
  double ta_i();
  double ta_v();
  double te_i();
  double te_v();
};