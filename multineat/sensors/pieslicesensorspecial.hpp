#ifndef HEADER_PIESLICESENSORSPECIAL
#define HEADER_PIESLICESENSORSPECIAL

// Second PieSliceSensor class to detect "special" enemies that require different behavior than usual
// For example jumping spikeys, can't be killed, shouldn't be jumped over
// Simply exist so the networks can react differently to these enemies when encountered

#include "pieslicesensor.hpp"
#include "../parameters.hpp"

class PieSliceSensorSpecial : public PieSliceSensor {
 friend class ExperimentParameterParser;
  
public:  
  PieSliceSensorSpecial(double angle1, double angle2) : PieSliceSensor(angle1, angle2) 
  {
    offset = rotate_point(Vector(ExperimentParameters::RADIUS_PIESLICE_SENSORS_SPECIAL, 0), angle1);
    offset2 = rotate_point(Vector(ExperimentParameters::RADIUS_PIESLICE_SENSORS_SPECIAL, 0), angle2);
  }
};

#endif /* HEADER_PIESLICESENSORSPECIAL */
/* EOF */