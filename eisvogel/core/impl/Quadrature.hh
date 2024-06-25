#pragma once

namespace Quadrature {

  struct TrapezoidalRule {

    // Calculates the weight of each sample for a trapezoidal quadrature in the
    // range [`ind_a`, `ind_b`) and `num_samples` is the total number of
    // quadrature points
    template <class IteratorT>
    static void fill_weights(std::size_t ind_a, std::size_t ind_b,
                             std::size_t num_samples, IteratorT weights_begin,
                             IteratorT weights_end)
    {

      // Make sure there is enough space
      assert((std::size_t)(weights_end - weights_begin) == ind_b - ind_a);

      // quadrature weights are 1.0 by default ...
      std::fill(std::execution::unseq, weights_begin, weights_end, 1.0);

      // ... except for the sample at the very beginning ...
      if (ind_a == 0) {
        [[unlikely]];
        *weights_begin = 0.5;
      }

      // ... or at the very end ...
      if (ind_b == num_samples) {
        [[unlikely]] * (weights_end - 1) = 0.5;
      }
    }
  };
} // namespace Quadrature
