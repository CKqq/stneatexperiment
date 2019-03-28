#ifndef PARAMETERS
#define PARAMETERS

#include <string>

using namespace std;

class ExperimentParameters {
public:
  static string path_to_level;
  
  static int num_cores;
  static int max_gens;
  
  static bool using_seed;
  static int seed;
  
  static string pop_filename;
  static string param_filename;
  static string experimentparam_filename;
  
  static int autosave_interval;
  static bool autosave_best;
  
  static int num_hidden_start_neurons;
  
  static int num_hidden_start_neurons_cppn;
  
  static bool novelty_search;
  static bool hyperneat;
  
  static int TILE_WIDTH;
  
  static int AMOUNT_RANGE_SENSORS;
  
  static int AMOUNT_DEPTH_SENSORS;
  static int SPACING_DEPTH_SENSORS;
  
  // This should be an odd number so we get a center/front pie slice
  static int AMOUNT_PIESLICE_SENSORS;
  static int RADIUS_PIESLICE_SENSORS;
  
  static int RADIUS_PIESLICE_SENSORS_SPECIAL;
};

#endif