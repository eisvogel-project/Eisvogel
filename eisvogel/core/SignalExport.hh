#pragma once

#include <vector>
#include <string>
#include "Common.hh"

void ExportSignal(std::vector<scalar_t> signal_times,
                  std::vector<scalar_t> signal_values, std::string outpath,
                  int time_prec = 3, int value_prec = 15);

#include "SignalExport.hxx"
