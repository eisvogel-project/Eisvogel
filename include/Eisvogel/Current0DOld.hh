#ifndef __CURRENT_0D_HH
#define __CURRENT_0D_HH

#include <vector>
#include "CoordUtils.hh"

#pragma message("Current0DOld is deprecated, please move to the new code.")

class Current0D {

public:

  Current0D() : m_points({}), m_charges({}) { };

  /**
   * Create a new particle trajectory as a collection of straight-line segments, with a certain
   * amount of charge flowing along each segment.
   *
   * @param points   Vector of (t, x, y, z) coordinates defining the start- and endpoints of the individual straight-line segments.
   *                 A vector of length `N` defines `N-1` line segments.
   * @param charges  Vector of length `N-1` containing charges for each line segment.
   */  
  Current0D(std::vector<CoordVector>&& points, std::vector<scalar_t>&& charges) : m_points(points), 
										  m_charges(charges) { 
    if(m_points.size() != m_charges.size() + 1) {
      throw;
    }
  };

  std::size_t number_segments() const {return m_charges.size();}

  const scalar_t& GetCharge(std::size_t ind) const {return m_charges[ind];}
  const FieldVector& GetPoint(std::size_t ind) const {return m_points[ind];}

private:

  std::vector<CoordVector> m_points;
  std::vector<scalar_t> m_charges;

};

#endif
