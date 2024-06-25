LineCurrentSegment::LineCurrentSegment(XYZCoordVector&& start_pos,
                                       XYZCoordVector&& end_pos,
                                       scalar_t start_time, scalar_t end_time,
                                       scalar_t charge)
    : start_pos(start_pos)
    , end_pos(end_pos)
    , start_time(start_time)
    , end_time(end_time)
    , charge(charge)
{

  if (start_time >= end_time) {
    throw std::runtime_error(
        "Error: have line current segment that spans zero or negative time "
        "period!");
  }
}
