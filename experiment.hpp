#ifndef EXPERIMENT_HEADER
#define EXPERIMENT_HEADER

#include <string>
#include <vector>

#include "multineat/Population.h"
#include "multineat/Genome.h"
#include "multineat/Species.h"
#include "multineat/Parameters.h"

using namespace std;
using namespace NEAT;

class Experiment {
public: 
  static void parse_files(string filename);
  
public:
  static int max_db_retry;
  static int db_sleeptime;
  
  Experiment();
  ~Experiment();
    
  void run_evolution();
  static void parse_commandline_arguments(int argc, char **argv);
  static void parse_experiment_parameters();
  static void fix_filenames();
  static Parameters init_params();
  
  static Parameters params;
  
  static string path_to_supertux_executable;
  static string path_to_level;
  
  static int num_cores;
  static int max_gens;
  
  static bool using_seed;
  static int seed;
  
  static string pop_filename;
  static string param_filename;
  static string experimentparam_filename;
  static string temp_pop_filename;
  static string db_filename;
  
  static int autosave_interval;
  
  static int num_range_sensors;
  static int num_depth_sensors;
  static int num_pieslice_sensors;
  
  static int num_hidden_start_neurons;
  
  
  static int busy_handler(void* data, int retry);
  static int select_handler(void* data, int argc, char** argv, char** colNames);
  
  static int cur_gen;
  
  static int top_genome_id;
  static double top_fitness;
  
  
  double xortest(Genome& g);
  bool constraints(Genome& g);
  
private:  
  // Static for db handler access
  static vector<Genome*> genomes;

  void refresh_genome_list();
  
  void save_pop(string filename);
  
  void start_processes();
  void update_db();
  void set_fitness_values();
  
  Genome start_genome;
  Population pop;
  
  static void print_usage();
};

#endif