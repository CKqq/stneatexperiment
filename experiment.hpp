#ifndef EXPERIMENT_HEADER
#define EXPERIMENT_HEADER

#include <string>
#include <vector>

#include "multineat/Population.h"
#include "multineat/Genome.h"
#include "multineat/Species.h"
#include "multineat/Parameters.h"
#include "multineat/Behavior.h"
#include "multineat/sensors/sensormanager.hpp"

#include <cxxabi.h>
#include <execinfo.h>
#include <stdio.h>

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
  
  static int busy_handler(void* data, int retry);
  static int select_handler(void* data, int argc, char** argv, char** colNames);
  static int select_handler_ns(void* data, int argc, char** argv, char** colNames);
  
  static double sparseness(Genome* g);
  void update_sparsenesses();
  
  static string temp_pop_filename;
  static string db_filename;
  
  static int cur_gen;
  
  static bool new_top;
  static int top_genome_gen_id;
  static int top_genome_id;
  static double top_fitness;
  
  static double evaluation_time;

  
private:
  // Static for db handler access
  void generate_substrate();
  
  static vector<Genome*> genomes;
  
  // Archive for Novelty Search
  static vector<PhenotypeBehavior> archive;

  void refresh_genome_list();
  
  void save_pop(string filename);
  
  void start_processes();
  void init_db();
  void update_db();
  void set_fitness_values();
  void update_gen_info();
  
  Genome start_genome;
  Population pop;
  
  Substrate substrate;
  
  static int num_winner_genomes;
  
  static float fitness_sum;
  static float airtime_sum;
  static float groundtime_sum;
  static int jump_sum;
  
  static void print_usage();
};

#endif