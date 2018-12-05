#include "sensormanager.hpp"

#include <math.h>

SensorManager* SensorManager::instance;


SensorManager::SensorManager() : 
  cur_sensors(new std::vector<std::shared_ptr<Sensor>>),
  cur_rangefinder_sensors(new std::vector<std::shared_ptr<RangeFinderSensor>>),
  cur_depthfinder_sensors(new std::vector<std::shared_ptr<DepthFinderSensor>>),
  cur_pieslice_sensors(new std::vector<std::shared_ptr<PieSliceSensor>>)
{
}

SensorManager::~SensorManager()
{
  
}

void SensorManager::initSensors()
{
  clearSensors();
  create_rangefinder_sensors();
  create_depthfinder_sensors();
  create_pieslice_sensors();
}

void SensorManager::clearSensors()
{
  cur_sensors->clear();
  cur_rangefinder_sensors->clear();
  cur_depthfinder_sensors->clear();
  cur_pieslice_sensors->clear();
}

std::vector<std::shared_ptr<Sensor>>* SensorManager::get_cur_sensors()
{
  return cur_sensors;
}


int SensorManager::get_total_sensor_count()
{
  return ExperimentParameters::AMOUNT_RANGE_SENSORS + ExperimentParameters::AMOUNT_DEPTH_SENSORS + ExperimentParameters::AMOUNT_PIESLICE_SENSORS * 2;
}


void SensorManager::create_rangefinder_sensors()
{
  int spacing = 0;
  
  for (int i = 0; i < ExperimentParameters::AMOUNT_RANGE_SENSORS; i++) {
    std::shared_ptr<RangeFinderSensor> sensor = std::make_shared<RangeFinderSensor>(spacing);
    cur_sensors->push_back(sensor);
    cur_rangefinder_sensors->push_back(sensor);
    spacing -= ExperimentParameters::TILE_WIDTH;
  }
}

void SensorManager::create_depthfinder_sensors()
{
  int spacing = 0;
  
  for (int i = 0; i < ExperimentParameters::AMOUNT_DEPTH_SENSORS; i++) {
    std::shared_ptr<DepthFinderSensor> sensor = std::make_shared<DepthFinderSensor>(spacing);
    cur_sensors->push_back(sensor);
    cur_depthfinder_sensors->push_back(sensor);
    spacing += ExperimentParameters::SPACING_DEPTH_SENSORS;
  }
}

void SensorManager::create_pieslice_sensors()
{
  double cur_angle = - M_PI / 2;
  double angle_step = M_PI / ExperimentParameters::AMOUNT_PIESLICE_SENSORS;
  
  for (int i = 0; i < ExperimentParameters::AMOUNT_PIESLICE_SENSORS; i++) {
    std::shared_ptr<PieSliceSensor> sensor = std::make_shared<PieSliceSensor>(cur_angle, cur_angle + angle_step);
    cur_sensors->push_back(sensor);
    cur_pieslice_sensors->push_back(sensor);
    
    sensor = std::make_shared<PieSliceSensorSpecial>(cur_angle, cur_angle + angle_step);
    cur_sensors->push_back(sensor);
    cur_pieslice_sensors->push_back(sensor);
    
    cur_angle += angle_step;
  }
}
