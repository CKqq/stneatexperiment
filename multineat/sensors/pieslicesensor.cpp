#include "pieslicesensor.hpp"

#include <math.h>

PieSliceSensor::PieSliceSensor(double angle1, double angle2) : Sensor(0, 0), 
angle1(angle1 < 0 ? angle1 + 2 * M_PI : angle1),
angle2(angle2 < 0 ? angle2 + 2 * M_PI : angle2)
{
  offset = rotate_point(Vector(SensorManager::RADIUS_PIESLICE_SENSORS, 0), angle1);
  offset2 = rotate_point(Vector(SensorManager::RADIUS_PIESLICE_SENSORS, 0), angle2);
}

Vector PieSliceSensor::rotate_point(const Vector& point, double rad)
{
  return Vector(point.x * cos(rad) - point.y * sin(rad), point.y * cos(rad) + point.x * sin(rad));
}

double PieSliceSensor::get_rad(const Vector& v1, const Vector& v2)
{
  double angle = atan2(v1.x * v2.y + v1.y * v2.x, v1.x * v2.x - v1.y * v2.y);
  
  if (angle < 0)
    angle += 2 * M_PI;
  
  return angle;
}

double PieSliceSensor::get_dist(const Vector& v1, const Vector& v2)
{
  return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y));
}

