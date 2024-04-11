#pragma once

#include "Vector.hh"

struct LineCurrentSegment {
  
  LineCurrentSegment(XYZCoordVector&& start_pos, XYZCoordVector&& end_pos, scalar_t start_time, scalar_t end_time, scalar_t charge);

  XYZCoordVector start_pos;
  XYZCoordVector end_pos;
  scalar_t start_time;
  scalar_t end_time;
  scalar_t charge;  
};
