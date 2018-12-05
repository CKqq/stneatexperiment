#ifndef HEADER_PIESLICESENSOR
#define HEADER_PIESLICESENSOR

#include "sensor.hpp"

class PieSliceSensor : public Sensor {
  
protected:  
  double angle1;
  double angle2;
  
  // With pie slice sensors, the offset values are only used for drawing
  // Also we need a second one since a pie slice is represented by two lines
  
  Vector offset2;
  
public:
  PieSliceSensor(double angle1, double angle2);
  
  Vector get_offset2() { return offset2; }
  
  // Rotates a point by given rad
  static Vector rotate_point(const Vector& point, double rad);
  
  // Returns angle between two points (or vectors) with the x axis as base
  static double get_rad(const Vector& v1, const Vector& v2);
  
  static double get_dist(const Vector& v1, const Vector& v2);
};

#endif /* HEADER_DEPTHFINDERSENSOR */
/* EOF */