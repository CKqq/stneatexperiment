#ifndef HEADER_RANGEFINDERSENSOR
#define HEADER_RANGEFINDERSENSOR

#include "sensor.hpp"

class RangeFinderSensor : public Sensor {
  
protected:
  static int length;
public:
  RangeFinderSensor(int offsetY);
};

#endif /* HEADER_RANGEFINDERSENSOR */
/* EOF */