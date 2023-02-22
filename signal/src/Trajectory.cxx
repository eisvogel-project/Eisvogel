#include "Trajectory.hh"

void Trajectory::AddPoint(CoordVector& point) {
  m_points.push_back(point);
}

void Trajectory::AddPoint(CoordVector&& point) {
  m_points.push_back(std::move(point));
}

CoordVector& Trajectory::operator()(std::size_t ind) {
  return m_points.at(ind);
}

const CoordVector& Trajectory::operator()(std::size_t ind) const {
  return m_points.at(ind);
}
