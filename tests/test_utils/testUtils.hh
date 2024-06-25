#include <cmath>
#include <format>

namespace TestUtils {

  template <typename T>
  bool close_match(T a, T b, T rel_th, T abs_th)
  {

    if ((a != 0.0) && (b != 0.0)) {
      return std::fabs((a - b) / b) < rel_th;
    }

    if (a == 0.0) {
      return std::fabs(b) < abs_th;
    }

    if (b == 0.0) {
      return std::fabs(a) < abs_th;
    }

    throw std::logic_error("Error: should never encounter this.");
  }

  // Checks ||sig_a - sig_b|| / min(||sig_a||, ||sig_b||) < rel_th
  template <typename T>
  bool signals_close_match(std::vector<T>& sig_a, std::vector<T>& sig_b,
                           T rel_th, bool verbose = true)
  {
    assert(sig_a.size() == sig_b.size());

    // Energies of the two input signals as well as their difference signal
    T sig_a_en = 0.0;
    T sig_b_en = 0.0;
    T sig_diff_en = 0.0;

    if (verbose) {
      std::cout << "---------------------------------------------------"
                << std::endl;
    }

    for (std::size_t i = 0; i < sig_a.size(); i++) {
      sig_a_en += std::pow(sig_a[i], 2.0);
      sig_b_en += std::pow(sig_b[i], 2.0);
      sig_diff_en += std::pow(sig_a[i] - sig_b[i], 2.0);

      if (verbose) {
        std::cout << std::format("A[{}] = {}, B[{}] = {}, (A-B)[{}] = {}", i,
                                 sig_a[i], i, sig_b[i], i, sig_a[i] - sig_b[i])
                  << std::endl;
      }
    }

    if (verbose) {
      std::cout << std::format("A = {}, B = {}, A - B = {}", sig_a_en, sig_b_en,
                               sig_diff_en)
                << std::endl;
      std::cout << "---------------------------------------------------"
                << std::endl;
    }

    // Test condition compares the energy of the difference signal to the
    // smaller of the two input energies
    return sig_diff_en / std::min(sig_a_en, sig_b_en) < rel_th;
  }
} // namespace TestUtils
