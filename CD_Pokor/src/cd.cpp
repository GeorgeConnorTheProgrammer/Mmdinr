#include <iostream>
#include <iomanip>
#include <numeric>
#include <array>

#define _USE_MATH_DEFINES
#include <cmath>

#include "cd.hpp"

//-----------------------------------------------------------------
// Simulation Functions
//-----------------------------------------------------------------

void CDState::Init()
{

}

void CDState::Step(double dt)
{
  double dCi1 = 0.0 * dt; // TODO - Eq 3a in Pokor
  double dCv1 = 0.0 * dt; // ^

  i_concentrations_swap[0] += dCi1;
  v_concentrations_swap[0] += dCv1;

  for (int n = 1; n <= NUM_CLUSTER_SIZES; ++n) // Compute the change in concentration for each cluster size
  {
    i_concentrations_swap[n] += dCi(n, dt);
    v_concentrations_swap[n] += dCv(n, dt);
  }

  rho += dRho(dt);

  std::swap(i_concentrations, i_concentrations_swap); // TODO - These two lines do not need to be O(NUM_CLUSTER_SIZES)
  std::swap(v_concentrations, v_concentrations_swap);
}

double CDState::C_i(int n)
{
  if (n < 1 || n > NUM_CLUSTER_SIZES) return 0.0;
  return i_concentrations[n - 1];
}

double CDState::C_v(int n)
{
  if (n < 1 || n > NUM_CLUSTER_SIZES) return 0.0;
  return v_concentrations[n - 1];
}

double CDState::dCi1(double dt)
{
  double D_i = D_0i * std::exp(-E_mi / (k * T));
  double D_v = D_0v * std::exp(-E_mv / (k * T));

  double R_iv = 4 * M_PI * (D_i + D_v) * r_iv; // TODO - Pokor Eq 3d

  return (G_i(1) - R_iv * C_i(1) * C_v(1) - C_i(1) / ta_gbi() - C_i(1) / ta_i() - C_i(1) / ta_i() + 1 / te_i()) * dt;
}

double CDState::dCv1(double dt)
{
  double D_i = D_0i * std::exp(-E_mi / (k * T));
  double D_v = D_0v * std::exp(-E_mv / (k * T));

  double R_iv = 4 * M_PI * (D_i + D_v) * r_iv;
  return (G_i(1) - R_iv * C_i(0) * C_v(0)) * dt;
}

double CDState::dCi(int n, double dt)
{
  return (G_i(n + 1) + a_i(n + 1) * C_i(n - 1) - b_i(n) * C_i(n) + c_i(n) * C_i(n + 1)) * dt;
}

double CDState::dCv(int n, double dt)
{
  return (G_v(n + 1) + a_v(n + 1) * C_v(n - 1) - b_v(n) * C_v(n) + c_v(n) * C_v(n + 1)) * dt;
}

double CDState::dRho(double dt) {
  return 1.0; // TODO - Sakaguchi Eq 3.14
}

double CDState::G_i(int n)
{
  // See Table 5 in Pokor
  const double nu = 0.3;
  const double G_dpa = 2.9 * std::pow(10, -7);
  const double fi2 = 0.5;
  const double fi3 = 0.2;
  const double fi4 = 0.06;

  switch(n)
  {
    case 1:
      return nu * G_dpa * (1 - fi2 - fi3 - fi4);
    case 2:
      return nu * G_dpa * fi2;
    case 3:
      return nu * G_dpa * fi3;
    case 4:
      return nu * G_dpa * fi4;
    default:
      return 0;
  }
}

double CDState::G_v(int n)
{
  // See Table 5 in Pokor
  const double nu = 0.3;
  const double G_dpa = 2.9 * std::pow(10, -7);
  const double fv2 = 0.06;
  const double fv3 = 0.03;
  const double fv4 = 0.02;

  switch(n)
  {
    case 1:
      return nu * G_dpa * (1 - fv2 - fv3 - fv4);
    case 2:
      return nu * G_dpa * fv2;
    case 3:
      return nu * G_dpa * fv3;
    case 4:
      return nu * G_dpa * fv4;
    default:
      return 0;
  }
}

double CDState::a_i(int n)
{
  return beta_iv(n + 1) * C_v(1) + alpha_ii(n + 1);
}

double CDState::b_i(int n)
{
  return beta_iv(n + 1) * C_v(1) + beta_ii(n) * C_i(1) + alpha_ii(n);
}

double CDState::c_i(int n)
{
  return beta_ii(n - 1) * C_i(1);
}

double CDState::a_v(int n)
{
  return beta_vi(n + 1) * C_i(1) + alpha_vv(n + 1);
}

double CDState::b_v(int n)
{
  return beta_vi(n + 1) * C_i(1) + beta_vv(n) * C_v(1) + alpha_vv(n);
}

double CDState::c_v(int n)
{
  return beta_vv(n - 1) * C_v(1);
}

double CDState::alpha_ii(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}
double CDState::alpha_vv(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}

double CDState::beta_iv(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}
double CDState::beta_vi(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}
double CDState::beta_ii(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}
double CDState::beta_vv(int n) {
  return 1.0; // TODO - Pokor Eq 4a-f
}

double CDState::ta_gbi()
{
  return 1.0; // TODO - Pokor Eq 3f
}

double CDState::ta_gbv()
{
  return 1.0; // TODO - Pokor Eq 3f
}

double CDState::ta_i()
{
  return 1.0; // TODO - Pokor Eq 3c
}

double CDState::ta_v()
{
  return 1.0; // TODO - Pokor Eq 3c
}

double CDState::te_i()
{
  return 1.0; // TODO - Pokor Eq 3b
}

double CDState::te_v()
{
  return 1.0 // TODO - Pokor Eq 3b
}