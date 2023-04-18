#ifndef __CURRENT_0D_HH
#define __CURRENT_0D_HH

#include <vector>
#include "CoordUtils.hh"

class Current0D {

public:

  Current0D() : m_points({}), m_charges({}) { };
  Current0D(std::vector<CoordVector>&& points, std::vector<scalar_t>&& charges) : m_points(points), 
										  m_charges(charges) { 
    if(m_points.size() != m_charges.size() + 1) {
      throw;
    }
  };

  void AddPoint(CoordVector& point, FieldVector& current);
  void AddPoint(CoordVector&& point, FieldVector&& current);

  std::size_t number_segments() const {return m_charges.size();}

  const scalar_t& GetCharge(std::size_t ind) const {return m_charges[ind];}
  const FieldVector& GetPoint(std::size_t ind) const {return m_points[ind];}

private:

  std::vector<CoordVector> m_points;
  std::vector<scalar_t> m_charges;

};

#endif
