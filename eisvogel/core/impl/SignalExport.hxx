#include "SignalExport.hh"
#include <fstream>
#include <iomanip>

void ExportSignal(std::vector<scalar_t> signal_times,
                  std::vector<scalar_t> signal_values, std::string outpath,
                  int time_prec, int value_prec)
{

  std::ofstream outfile;
  outfile.open(outpath);
  for (int ind = 0; ind < signal_times.size(); ind++) {
    outfile << std::fixed << std::showpoint << std::setprecision(time_prec);
    outfile << signal_times[ind];
    outfile << ", ";
    outfile << std::scientific << std::showpoint
            << std::setprecision(value_prec);
    outfile << signal_values[ind] << std::endl;
  }
  outfile.close();
}
