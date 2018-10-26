#ifndef HEADER_SENSOR
#define HEADER_SENSOR

#include <vector>

class Vector{
public:
  double x;
  double y;
  
  Vector() : x(0), y(0)
  {
  }
  
  Vector(double x, double y) : x(x), y(y)
  {
  }
};

class Sensor{
  
public:
  Sensor(int offsetX, int offsetY);
  
protected:
  Vector offset;
  double value;
  
public:
  Vector get_offset() { return offset; };
  double getValue() { return value; };
};

#endif /* HEADER_SENSOR */
/* EOF */