#include <iostream>
#include <fstream>
#include <numeric>
#include <array>
#include "TGraph.h"
#include "TCanvas.h"
#include "TGButton.h"
#include "TGClient.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TF1.h"

#define _USE_MATH_DEFINES
#include<cmath>

#include "../vendor/nlohmann/json.hpp"

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
  void run_model(double sample_interval, double dt, double endtime, double Temp, int K_0_exp, int C_s_exp, TGraph* g1, TGraph* g2) 
  {
		// running variable init
		// Calculating the D_i and D_v.
		double D_i = D_0i * std::exp(-(E_mi/(k*Temp))); // Interstitial diffusion coefficient
		double D_v = D_0v * std::exp(-(E_mv/(k*Temp))); // Vacancy diffusion coefficient		
		
		double K_0 = std::pow(10,K_0_exp); //defect production rate
		double C_s = std::pow(10,C_s_exp); //sink concentration

		double K_iv = 4.0 * M_PI * r_iv * (D_i + D_v); // vancancy-interstitial recombination rate coeff
		double K_is = 4.0 * M_PI * r_is * D_i; // interstitial-sink reaction rate coeff.
		double K_vs = 4.0 * M_PI * r_vs * D_v; // vancancy-sink reaction rate coeff.

    double sample_counter = sample_interval;
    for (double t = 0.0; t < endtime; t += dt) 
    {
      double dC_i = (K_0 - K_iv * C_i * C_v - K_is * C_i * C_s) * dt;
      double dC_v = (K_0 - K_iv * C_i * C_v - K_vs * C_v * C_s) * dt;

			if (std::isinf(dC_i) || std::isinf(dC_v) || std::isnan(dC_i) || std::isnan(dC_v))
			{
				std::cerr << "Limit reached stopping model.." << std::endl;
				return; 
			}

      C_i += dC_i;
      C_v += dC_v;

      sample_counter += dt;
      if (sample_counter >= sample_interval)
      {
				// double sink_diff = K_is * C_i * C_s - K_vs * C_v * C_s;
        sample_counter = 0.0;
				g1->AddPoint(t, C_i);
				g2->AddPoint(t, C_v);
      }
    }

    std::cerr << "Finished simulation" << std::endl;
  }

int main(int argc, char** argv)
{
  // Define the application start point and the name <- args of main.
  TApplication app("app", &argc, argv);
  
  // Simulation ..
  const std::string config_file_name = argc < 2 ? "config.json" : argv[1];

  std::ifstream config_file(config_file_name);
  if (!config_file.good()) 
  {
    std::cerr << "Could not open " << config_file_name << std::endl;
    return 1;
  }

  nlohmann::json config = nlohmann::json::parse(config_file);
  double time = config["total_time_seconds"];
  double dt = config["dt_seconds"];
  double Temp = config["temperature_kelvin"];
  int K_0_exp = config["K_0_exp"];
  int C_s_exp = config["C_s_exp"];

  double sample_interval = config["sample_interval"];
  // Define the 'canvas' aka main window
  TCanvas* c = new TCanvas("c", "Sim", 0, 0, 800, 600);
  
  // c->SetLogy();
  
  // TGTextButton *but1 = new TGTextButton("Run Simulation","",.05,.8,.45,.88);
  // but1->Draw();
  
  // Define the graphs
  TGraph* ci_graph = new TGraph();
  TGraph* cv_graph = new TGraph();

  ci_graph->SetMarkerColor(kRed);
  cv_graph->SetMarkerColor(kBlue);

	ci_graph->SetTitle("C_i");
	cv_graph->SetTitle("C_v");

  TMultiGraph* mg = new TMultiGraph();
  
  run_model(sample_interval, dt, time, Temp, K_0_exp, C_s_exp, ci_graph, cv_graph);
  
  mg->Add(ci_graph);
  mg->Add(cv_graph);

  mg->SetTitle("MFRT");
  mg->Draw("ALP");
  
  // Graph config / color
  c->Modified(); c->Update();
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();

  // Term app on window close
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  
  // Simulation ..

  // App run
  app.Run();


  return 0;
}
