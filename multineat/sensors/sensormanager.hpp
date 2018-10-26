#ifndef SENSOR_MANAGER
#define SENSOR_MANAGER

#include "sensor.hpp"
#include "rangefindersensor.hpp"
#include "depthfindersensor.hpp"
#include "pieslicesensor.hpp"

#include <memory>
#include <vector>
#include <iostream>

class PieSliceSensor;

class SensorManager {
  friend class Experiment;
  friend class PieSliceSensor; // Dirty Workaround, generate Offsets in this class
  
private:
  static int TILE_WIDTH;
  
  static int AMOUNT_RANGE_SENSORS;
  
  static int AMOUNT_DEPTH_SENSORS;
  static int SPACING_DEPTH_SENSORS;
  
  // This should be an odd number so we get a center/front pie slice
  static int AMOUNT_PIESLICE_SENSORS;
  static int RADIUS_PIESLICE_SENSORS;
public:
  SensorManager();
  ~SensorManager();
  
  void initSensors();
  void clearSensors();
  
  std::vector<std::shared_ptr<Sensor>>* get_cur_sensors();
  
  std::vector<std::shared_ptr<RangeFinderSensor>>* get_cur_rangefinder_sensors() { return cur_rangefinder_sensors; };
  std::vector<std::shared_ptr<DepthFinderSensor>>* get_cur_depthfinder_sensors() { return cur_depthfinder_sensors; };
  std::vector<std::shared_ptr<PieSliceSensor>>*	   get_cur_pieslice_sensors()    { return cur_pieslice_sensors; };
  
  static SensorManager* instance;
  
  static int get_total_sensor_count();
  
public:
  std::vector<std::shared_ptr<Sensor>>* cur_sensors;
  std::vector<std::shared_ptr<RangeFinderSensor>>* cur_rangefinder_sensors;
  std::vector<std::shared_ptr<DepthFinderSensor>>* cur_depthfinder_sensors;
  std::vector<std::shared_ptr<PieSliceSensor>>*	   cur_pieslice_sensors;
  void create_rangefinder_sensors();
  void create_depthfinder_sensors();
  void create_pieslice_sensors();
};

#endif /* SENSOR_MANAGER */
/* EOF */