#include "shower_1D.h"
#include <array>
#include <iostream>
#include "shower_creator.h"
#include "TCanvas.h"
#include "TGraph.h"
int main()
{
  std::array<float, 3> pos = {20, 20, 20};
  showers::ShowerCreator shower_creator(
      "/home/welling/RadioNeutrino/scripts/Eisvogel/shower_profile/"
      "shower_file");
  double pi = acos(-1);

  showers::Shower1D new_shower =
      shower_creator.create_shower(pos, 2.e16, .2 * pi, 0.3 * pi, 1);
  std::vector<double> shower_x;
  std::vector<double> shower_y;
  std::vector<double> shower_z;
  std::vector<double> shower_t;
  std::vector<double> current_x;
  std::vector<double> current_y;
  std::vector<double> current_z;

  new_shower.get_current(1.e-9, &shower_t, &shower_x, &shower_y, &shower_z,
                         &current_x, &current_y, &current_z);

  TCanvas* c = new TCanvas("c", "Something", 0, 0, 800, 600);
  auto graph1 = new TGraph(shower_t.size(), &shower_t[0], &current_z[0]);
  c->SetLogy();
  graph1->Draw();
  c->Print("plot_shower_profile.png");
}
