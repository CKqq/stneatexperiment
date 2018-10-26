#ifndef HEADER_DEPTHFINDERSENSOR
#define HEADER_DEPTHFINDERSENSOR

#include "sensor.hpp"

class DepthFinderSensor : public Sensor {
  friend class ExperimentParameterParser;
  
protected:
  static int length;
public:
  DepthFinderSensor(int offsetX);
  
};

#endif /* HEADER_DEPTHFINDERSENSOR */
/* EOF */