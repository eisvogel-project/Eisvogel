#ifndef __TRAJECTORY_HH
#define __TRAJECTORY_HH

#include <vector>
#include "CoordUtils.hh"

class Trajectory {

public:
  
  Trajectory() : m_points({}) { };
  Trajectory(std::vector<CoordVector>&& points) : m_points(points) { };
  
  void AddPoint(CoordVector& point);
  void AddPoint(CoordVector&& point);

  std::size_t size() const {return m_points.size();}

  CoordVector& operator()(std::size_t ind);
  const CoordVector& operator()(std::size_t ind) const;

  auto begin() {return m_points.begin();}
  auto begin() const {return m_points.cbegin();}
  auto cbegin() {return m_points.cbegin();}
  auto end() {return m_points.end();}
  auto end() const {return m_points.cend();}
  auto cend() {return m_points.cend();}

private:

  std::vector<CoordVector> m_points;

};

#endif
