#include "rangefindersensor.hpp"

int RangeFinderSensor::length = 256;

RangeFinderSensor::RangeFinderSensor(int offsetY) : Sensor(0, offsetY)
{
}