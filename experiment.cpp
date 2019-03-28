#include <iostream>
#include <string.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <chrono>

#include <sys/stat.h>

#include <sqlite3.h>

#include "experiment.hpp"

#include <signal.h>

string Experiment::temp_pop_filename = "temp_pop";

string Experiment::db_filename = "neat.db";

int Experiment::max_db_retry = 100;
int Experiment::db_sleeptime = 50;

string Experiment::path_to_supertux_executable = "./bin/build/supertux2";

int Experiment::cur_gen = 0;
bool Experiment::new_top = true;
int Experiment::top_genome_gen_id = 0;
int Experiment::top_genome_id = 0;
double Experiment::top_fitness = 0;

int Experiment::num_winner_genomes = 0;
float Experiment::fitness_sum = 0;
float Experiment::airtime_sum = 0;
float Experiment::groundtime_sum = 0;
int Experiment::jump_sum = 0;
double Experiment::evaluation_time = 0;


Parameters Experiment::params;

vector<Genome*> Experiment::genomes;
vector<PhenotypeBehavior> Experiment::archive;

int main(int argc, char **argv) {
  std::cout << "SuperTux + NEAT interface and experiment code by Christoph Kuhfuss 2018" << std::endl;
  std::cout << "Using the original SuperTux source and the MultiNEAT framework by Peter Chervenski (https://github.com/peter-ch/MultiNEAT)" << std::endl;
  
  Experiment::parse_commandline_arguments(argc, argv);
  Experiment::fix_filenames();
  Experiment::parse_experiment_parameters();
  Experiment::params = Experiment::init_params();
  
  Experiment experiment;
  
  experiment.run_evolution();
  
  return 0;
}


Experiment::Experiment() : 
start_genome(0, ExperimentParameters::AMOUNT_RANGE_SENSORS + ExperimentParameters::AMOUNT_DEPTH_SENSORS + ExperimentParameters::AMOUNT_PIESLICE_SENSORS * 2 + 1, ExperimentParameters::num_hidden_start_neurons, 
  6, false, UNSIGNED_SIGMOID, UNSIGNED_SIGMOID, 1, params),
  pop(strcmp(ExperimentParameters::pop_filename.c_str(), "") ? 
  Population(ExperimentParameters::pop_filename.c_str()) : Population(start_genome, params, true, 2.0, (ExperimentParameters::using_seed ? ExperimentParameters::seed : (int) time(0)))),
  total_gen_time(0)
{
  if (ExperimentParameters::hyperneat) 
  {
    generate_substrate();
  }
}

Experiment::~Experiment()
{
}

void Experiment::run_evolution() {
  std::remove(db_filename.c_str());
  
  init_db();
  
  for (int i = 1; i <= ExperimentParameters::max_gens; i++) {      
    std::cout << "Working on generation #" << i << "..." << std::endl;
    
    Experiment::new_top = false;
    
    num_winner_genomes = 0;
    
    fitness_sum = 0;
    airtime_sum = 0;
    groundtime_sum = 0;
    jump_sum = 0;
    
    evaluation_time = 0;
    
    cur_gen = i;
    refresh_genome_list();
    
    if (ExperimentParameters::novelty_search) {
      std::vector<PhenotypeBehavior>* behaviors = new std::vector<PhenotypeBehavior>();
      
      for (int i = 0; i < genomes.size(); i++) 
	behaviors->push_back(Behavior());
      
      pop.InitPhenotypeBehaviorData(behaviors, &archive);
    }
    
    // Prepare db file...
    update_db();
    
    // Run supertux processes...
    start_processes();
    
    // ... and get the fitness values from the db file
    // Also update info about the generation
    set_fitness_values();
    update_gen_info();
    
    std::cout << "Done. Top fitness: " << top_fitness << " by individual #" << top_genome_id << ", generation #" << top_genome_gen_id << std::endl;
    std::cout << "Evaluation time: " << (int) (evaluation_time / 60) << "min" << (int) evaluation_time % 60 << "sec" << std::endl;
    std::cout << num_winner_genomes << " genome(s) out of " << genomes.size() << " finished the level (" << num_winner_genomes  / (double) genomes.size() * 100 << "%)" << std::endl;
    
        
    if (ExperimentParameters::autosave_best && Experiment::new_top) {
      std::ostringstream ss;
      ss << "./neat_gen" << cur_gen;
      save_pop(ss.str().c_str());
    }
    else if (!ExperimentParameters::autosave_best && i % ExperimentParameters::autosave_interval == 0) {
      std::ostringstream ss;
      ss << "./neat_gen" << cur_gen;
      save_pop(ss.str().c_str());
    }
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    pop.Epoch();
    std::chrono::time_point<std::chrono::high_resolution_clock> finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;
    
    total_gen_time += elapsed.count();
  }
}

void Experiment::start_processes()
{
  // Save pop to single file for all child processes to read

  pop.Save(temp_pop_filename.c_str());
  
  // Distribute genomes as evenly as possible
  int remaining_genome_count = genomes.size();
  int cores_left = ExperimentParameters::num_cores;
  
  std::vector<int> startindices;
  std::vector<int> endindices;
  
  int curindex = 0;
  
  for (int j = 0; j < ExperimentParameters::num_cores; j++) {
    startindices.push_back(++curindex - 1);
    
    int single_genome_count = remaining_genome_count / cores_left;
    
    if (remaining_genome_count % cores_left) 
      single_genome_count++;
              
    remaining_genome_count -= single_genome_count;
    cores_left--;
    
    curindex += single_genome_count - 1;
    endindices.push_back(curindex - 1);
  }
  
  // Build the GNU parallel command, include all parameter files and the levelfile
  std::stringstream ss;
  ss << "parallel --no-notice --xapply ";
  ss << path_to_supertux_executable;
  ss << " --neat";
  
//   if (Experiment::novelty_search) ss << " --noveltysearch";
  if (ExperimentParameters::hyperneat) ss << " --hyperneat";
  
  // All child processes read from the same population file
  ss << " --popfile " << boost::filesystem::canonical(temp_pop_filename);
  
  // Only add experiment and MultiNEAT param files if they're configured...
  if (strcmp(ExperimentParameters::experimentparam_filename.c_str(), "")) 	ss << " --experimentparamfile " << ExperimentParameters::experimentparam_filename;
  if (strcmp(ExperimentParameters::param_filename.c_str(), "")) 		ss << " --paramfile " << ExperimentParameters::param_filename;
  
  ss << " --dbfile " << boost::filesystem::canonical(db_filename);
  ss << " --curgen " << cur_gen;
  
  ss <<	" --fromtogenome ::: ";
  
  for (std::vector<int>::iterator it = startindices.begin(); it != startindices.end(); ++it) {
    ss << *it << " ";
  }
  
  ss << "::: ";
  
  for (std::vector<int>::iterator it = endindices.begin(); it != endindices.end(); ++it) {
    ss << *it << " ";
  }
  
  ss << "::: '";
  
  ss << ExperimentParameters::path_to_level;
  
  ss << "'";
  
//   std::cout << ss.str() << std::endl;
  
  // Ready to rumble, don't forget to measure the time
  std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
  std::system(ss.str().c_str());
  std::chrono::time_point<std::chrono::high_resolution_clock> finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = finish - start;
  
  evaluation_time = elapsed.count();
  // Remove temporary population file. If we have to, we'll make a new one
  std::remove(temp_pop_filename.c_str());
}

void Experiment::init_db()
{
  sqlite3* db;  
  sqlite3_open("neat.db", &db);
  
  sqlite3_busy_handler(db, busy_handler, (void*) nullptr);
    
  char* err;
  std::stringstream ss;
  
  ss << "CREATE TABLE RUN_INFO (id INT PRIMARY KEY NOT NULL, num_genomes INT, avg_fitness REAL, time REAL, avg_airtime REAL, avg_groundtime REAL, avg_jumps REAL, top_fitness REAL, amt_winners INT);";
  
  sqlite3_exec(db, ss.str().c_str(), 0, 0, &err);
  
  sqlite3_close(db);
}

void Experiment::update_db()
{
  sqlite3* db;  
  sqlite3_open("neat.db", &db);
  
  sqlite3_busy_handler(db, busy_handler, (void*) nullptr);
    
  char* err;
  std::stringstream ss;
  
  ss.str("");
  
  ss << "CREATE TABLE GEN" << cur_gen << "(id INT PRIMARY KEY NOT NULL, fitness REAL, airtime REAL, groundtime REAL, num_jumps INT, qLeft REAL, qRight REAL, qUp REAL, qDown REAL, qJump REAL, qAction REAL, won INT);";
    
  sqlite3_exec(db, ss.str().c_str(), 0, 0, &err);
  
  for(std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    ss.str("");
    
    // Negative fitness if not evaluated yet. If fitness values stay negative, NEAT might just break
    // Last values are for ns behavior
    ss << "INSERT INTO GEN" << cur_gen << " VALUES(" << (*it)->GetID() << ", -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);";
    
    sqlite3_exec(db, ss.str().c_str(), 0, 0, &err);
  }

  sqlite3_close(db);
}


void Experiment::set_fitness_values()
{
  sqlite3* db;  
  sqlite3_open("neat.db", &db);
  
  sqlite3_busy_handler(db, busy_handler, (void*) nullptr);
  
  char* err;
  
  std::stringstream ss;
  ss << "SELECT * FROM gen" << cur_gen << ";";
  
  sqlite3_exec(db, ss.str().c_str(), (ExperimentParameters::novelty_search ? select_handler_ns : select_handler), 0, &err);
  
  if (ExperimentParameters::novelty_search) {
    for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
      if ((*it)->m_PhenotypeBehavior->m_Data[0].size() == 0) {
	std::cout << "ERROR OCCURED WITH GENOME #" << (*it)->GetID() << std::endl;
      }
    }
    update_sparsenesses();
  }

  sqlite3_close(db);
}

void Experiment::update_gen_info()
{
  sqlite3* db;  
  sqlite3_open("neat.db", &db);
  
  sqlite3_busy_handler(db, busy_handler, (void*) nullptr);
    
  char* err;
  std::stringstream ss;

  ss << "INSERT INTO RUN_INFO VALUES(" << cur_gen << ", " << genomes.size() << ", " << fitness_sum / genomes.size() << ", " << evaluation_time << ", " << (airtime_sum / (double) genomes.size()) << ", " << (groundtime_sum / (double) genomes.size()) << ", " << ((double) jump_sum) / genomes.size() << ", " << top_fitness << ", " << num_winner_genomes << ");";
    
  sqlite3_exec(db, ss.str().c_str(), 0, 0, &err);

  sqlite3_close(db);
}

void Experiment::refresh_genome_list()
{
  genomes.clear();
  
  for (unsigned int i = 0; i < pop.m_Species.size(); i++) {
    Species* cur = &pop.m_Species[i];
//     std::cout << "Species #" << cur->ID() << " has " << cur->m_Individuals.size() << " individuals" << std::endl;
    for (unsigned int j = 0; j < cur->m_Individuals.size(); j++) {
      cur->m_Individuals[j].m_PhenotypeBehavior = new Behavior();
      genomes.push_back(&cur->m_Individuals[j]);
    }
    
    // Probably unneccessary
    sort(genomes.begin(), genomes.end(), [ ](const Genome* lhs, const Genome* rhs)
    {
      return lhs->GetID() < rhs->GetID();
    });
  }
}

double Experiment::sparseness(Genome* g)
{
  std::vector<double> distances;
  
  double sum = 0;
  int amt = 0;
  
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    //std::cout << "Calculating distance between genomes #" << g->GetID() << " and #" << (*it)->GetID() << std::endl;
    amt++;
    distances.push_back(g->m_PhenotypeBehavior->Distance_To((*it)->m_PhenotypeBehavior));
  }
  
  for (std::vector<PhenotypeBehavior>::iterator it = archive.begin(); it != archive.end(); ++it) {
    amt++;
    distances.push_back(g->m_PhenotypeBehavior->Distance_To(&(*it)));
  }
  
  sort(distances.begin(), distances.end(), [ ](const double lhs, const double rhs)
  {
      return lhs < rhs;
  });
  
  for (int i = 0; i < params.NoveltySearch_K; i++) {
    sum += distances[i];
  }
  
//   std::cout << "Sparseness for genome #" << g->GetID() << " is " << sum / amt << std::endl;
 
  if (amt > 0) {
    double sparseness = sum / amt;
    
    if (sparseness > params.NoveltySearch_P_min) {
      archive.push_back(*(g->m_PhenotypeBehavior));
    }
    
    return sum / amt;
  }
  else
    return 0;
}

void Experiment::update_sparsenesses()
{
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it)
  {        
    (*it)->SetFitness(pop.ComputeSparseness(*(*it)));
    (*it)->SetEvaluated();
  }
}

void Experiment::generate_substrate()
{
  SensorManager sm;
  sm.initSensors();
  
  std::vector<std::vector<double>> input_coords;
  std::vector<std::vector<double>> hidden_coords;	// TODO: Add support for hidden neurons (third dimension, arrange in a circle)
  std::vector<std::vector<double>> output_coords;
  
  std::vector<std::shared_ptr<RangeFinderSensor>>* rfsensors = sm.get_cur_rangefinder_sensors();
  std::vector<std::shared_ptr<DepthFinderSensor>>* dfsensors = sm.get_cur_depthfinder_sensors();
  std::vector<std::shared_ptr<PieSliceSensor>>*    pssensors = sm.get_cur_pieslice_sensors();
  
  double inputZ = -1;
  double outputZ = 1;
  double hiddenZ = 0;
  
  // First, calculate maximum x and y coordinates so we can normalize substrate coordinates to (-1, 1)
  int maxX = 0;
  int maxY = 0;
  
  if (ExperimentParameters::num_hidden_start_neurons > 0) {
    std::vector<double> coords;
    if (ExperimentParameters::num_hidden_start_neurons == 1) {
      // If there is only one hidden start neuron, place it in the center and be done with it
      coords.push_back(0);
      coords.push_back(0);
      coords.push_back(hiddenZ);
      
      hidden_coords.push_back(coords);
    } else {
      // If there are more than one, arrange in a circle
      coords.push_back(1);
      coords.push_back(0);
      coords.push_back(hiddenZ);
      
      double cur_angle = 0;
      
      // How many neurons per full rotation?
      double angle_step = 2 * M_PI / ExperimentParameters::num_hidden_start_neurons;
      
      for (int i = 0; i < ExperimentParameters::num_hidden_start_neurons - 1; i++) {
	// Push back coordinates, rotate by the angle step each iteration
	hidden_coords.push_back(coords);
	
	std::vector<double> coords_new;
	
	coords_new.push_back(coords.at(0) * cos(angle_step) + coords.at(1) * sin(angle_step));
	coords_new.push_back(coords.at(1) * cos(angle_step) - coords.at(0) * sin(angle_step));
	coords_new.push_back(hiddenZ);

	coords = coords_new;
      }
    }
  }
  
  for (std::vector<std::shared_ptr<RangeFinderSensor>>::iterator it = rfsensors->begin(); it != rfsensors->end(); ++it)
  {
    Vector v = (*it)->get_offset();
    
    // Get middle point of the RangeFinderSensor
    int x = abs(v.x / 2.0);
    int y = abs(v.y);
    
    if (x > maxX)
      maxX = x;
    
    if (y > maxY) 
      maxY = y;
    
  }
  
  for (std::vector<std::shared_ptr<DepthFinderSensor>>::iterator it = dfsensors->begin(); it != dfsensors->end(); ++it)
  {
    Vector v = (*it)->get_offset();
    
    // Get middle point of the DepthFinderSensor
    int x = abs(v.x);
    int y = abs(v.y / 2.0);
    
    if (x > maxX) 
      maxX = x;
    
    if (y > maxY)
      maxY = y;
    
  }
  
  for (std::vector<std::shared_ptr<PieSliceSensor>>::iterator it = pssensors->begin(); it != pssensors->end(); ++it)
  {
    Vector v1 = (*it)->get_offset();
    Vector v2 = (*it)->get_offset2();
    
    // Get middle point of the PieSliceSensor
    // Divide by three because both vectors form a triangle with (0, 0) as the third point
    int x = abs((v1.x + v2.x) / 3.0);
    int y = abs((v1.y + v2.y) / 3.0);
    
    if (x > maxX)
      maxX = x;
    
    if (y > maxY)
      maxY = y;
    
  }
    
  for (std::vector<std::shared_ptr<RangeFinderSensor>>::iterator it = rfsensors->begin(); it != rfsensors->end(); ++it)
  {
    std::vector<double> coords;
    
    Vector v = (*it)->get_offset();
    
    // Normalize between (0, 1) and then between (-1, 1)
    coords.push_back(v.x / (double) maxX);
    coords.push_back(v.y / (double) maxY);
    
    coords.push_back(v.x);
    coords.push_back(v.y);
    
    if (ExperimentParameters::num_hidden_start_neurons > 0) coords.push_back(inputZ);
    
    input_coords.push_back(coords);
  }
  
  for (std::vector<std::shared_ptr<DepthFinderSensor>>::iterator it = dfsensors->begin(); it != dfsensors->end(); ++it)
  {
    std::vector<double> coords;
    
    Vector v = (*it)->get_offset();
    
    coords.push_back((v.x / (double) maxX));
    coords.push_back((v.y / (double) maxY));
    
    coords.push_back(v.x);
    coords.push_back(v.y);
    
    if (ExperimentParameters::num_hidden_start_neurons > 0) coords.push_back(inputZ);
    
    input_coords.push_back(coords);
  }
  
  for (std::vector<std::shared_ptr<PieSliceSensor>>::iterator it = pssensors->begin(); it != pssensors->end(); ++it)
  {
    std::vector<double> coords;
    
    Vector v1 = (*it)->get_offset();
    Vector v2 = (*it)->get_offset2();
    
    // Same procedure as above, third point is (0, 0)    
    coords.push_back(((v1.x + v2.x) / 3.0) / (double) maxX);
    coords.push_back(((v1.y + v2.y) / 3.0) / (double) maxY);
    if (ExperimentParameters::num_hidden_start_neurons > 0) coords.push_back(inputZ);
    
    input_coords.push_back(coords);
  }
  
  std::vector<double> bias;
  bias.push_back(0);
  bias.push_back(0);
  if (ExperimentParameters::num_hidden_start_neurons > 0) bias.push_back(inputZ);
  
  input_coords.push_back(bias);
    
  std::vector<double> output;
  
  // Arrange output substrate coordinates around Tux
  // UP should be above Tux, DOWN below, etc. with Tux "being" at (0, 0)
  
  // UP
  output.push_back(0);
  output.push_back(1);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);
  
  output_coords.push_back(output);
  output.clear();
  
  // DOWN
  output.push_back(0);
  output.push_back(-1);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);

  
  output_coords.push_back(output);
  output.clear();
  
  // LEFT
  output.push_back(-1);
  output.push_back(0);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);
  
  output_coords.push_back(output);
  output.clear();
  
  // RIGHT
  output.push_back(1);
  output.push_back(0);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);
  
  output_coords.push_back(output);
  output.clear();
  
  // JUMP
  output.push_back(0);
  output.push_back(0.8);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);
  
  output_coords.push_back(output);
  output.clear();
  
  // ACTION
  output.push_back(0);
  output.push_back(0.1);
  if (ExperimentParameters::num_hidden_start_neurons > 0) output.push_back(outputZ);
  
  output_coords.push_back(output);
  
  // Parameters taken from MultiNEAT's HyperNEAT example
  
  substrate = Substrate(input_coords, hidden_coords, output_coords);
  
  substrate.m_allow_hidden_hidden_links = false;
  substrate.m_allow_hidden_output_links = true;
  substrate.m_allow_input_hidden_links = true;
  substrate.m_allow_input_output_links = true;
  substrate.m_allow_looped_hidden_links = false;
  substrate.m_allow_looped_output_links = false;
  substrate.m_allow_output_hidden_links = false;
  substrate.m_allow_output_output_links = false;
  substrate.m_query_weights_only = true;
  substrate.m_max_weight_and_bias = 8;
  
  substrate.m_hidden_nodes_activation = ActivationFunction::SIGNED_SIGMOID;
  substrate.m_output_nodes_activation = ActivationFunction::UNSIGNED_SIGMOID;
  
  start_genome = Genome(0, substrate.GetMinCPPNInputs(), ExperimentParameters::num_hidden_start_neurons_cppn, 
	       substrate.GetMinCPPNOutputs(), false, TANH, TANH, 0, params);
  pop = strcmp(ExperimentParameters::pop_filename.c_str(), "") ? Population(ExperimentParameters::pop_filename.c_str()) : Population(start_genome, params, true, 1.0, (ExperimentParameters::using_seed ? ExperimentParameters::seed : (int) time(0)));
}


int Experiment::select_handler(void* data, int argc, char** argv, char** colNames)
{
  int id = std::stoi(argv[0]);
  double fitness = std::stod(argv[1]);
  
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    if ((*it)->GetID() == id) {
      (*it)->SetFitness((ExperimentParameters::novelty_search) ? sparseness(*it) : fitness);
      (*it)->SetEvaluated();
      if (fitness > top_fitness) {
	new_top = true;
	
	top_fitness = fitness;
	top_genome_id = id;
	top_genome_gen_id = cur_gen;
      }
      
      break;
    }
  }
  
  fitness_sum += std::stod(argv[1]);
  airtime_sum += std::stod(argv[2]);
  groundtime_sum += std::stod(argv[3]);
  jump_sum += std::stoi(argv[4]);
  
  if (std::stod(argv[11]) > 0) num_winner_genomes++;
  
  return 0;
}

int Experiment::select_handler_ns(void* data, int argc, char** argv, char** colNames)
{
  int id = std::stoi(argv[0]);

  double fitness = std::stod(argv[1]);
  
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    if ((*it)->GetID() == id) {
            
      for (int i = 5; i < 11; i++) {
	(*it)->m_PhenotypeBehavior->m_Data[0][i - 5] = std::stod(argv[i]);
      }
      
      if (fitness > top_fitness) {
	new_top = true;
	top_fitness = fitness;
	top_genome_id = id;
	top_genome_gen_id = cur_gen;
      }
      
      break;
    }
  }

  fitness_sum += std::stod(argv[1]);
  airtime_sum += std::stod(argv[2]);
  groundtime_sum += std::stod(argv[3]);
  jump_sum += std::stoi(argv[4]);
  
  if (std::stod(argv[11]) > 0) num_winner_genomes++;
  
  return 0;
}

void Experiment::parse_commandline_arguments(int argc, char** argv)
{
  for (int i = 0; i < argc; i++) {
    string arg = argv[i];
    
    if (arg == "--help")
    {
      print_usage();
      exit(0);
    }
    else if (arg == "--experimentparamfile")
    {
      if (i + 1 < argc)
      {
	ExperimentParameters::experimentparam_filename = argv[++i];
      }
      else
      {
	std::cout << "Please specify path to experiment parameter file after using --experimentparamfile" << std::endl;
      }
    }
    else if (arg == "--paramfile")
    {
      if (i + 1 < argc)
      {
	ExperimentParameters::param_filename = argv[++i];
      }
      else
      {
	std::cout << "Please specify path to MultiNEAT parameter file after using --paramfile" << std::endl;
      }
    }
    else if (arg == "--popfile")
    {
      if (i + 1 < argc)
      {
	ExperimentParameters::pop_filename = argv[++i];
      }
      else
      {
	std::cout << "Please specify path to population file after using --popfile" << std::endl;
      }
    }
    
    // We'll take seeds as a console argument for convenience,
    // i.e. start multiple experiments via shellscript
    else if (arg == "--seed")
    {
      int seed = 0;
      
      if (i + 1 < argc && sscanf(argv[i + 1], "%d", &seed) == 1) {
	ExperimentParameters::using_seed = true;
      } else {
	std::cout << "Please specify a proper seed after '--seed'" << std::endl;
      }
    }
    else if (arg == "--noveltysearch")
    {
      ExperimentParameters::novelty_search = true;
    }
    else if (arg == "--hyperneat")
    {
      ExperimentParameters::hyperneat = true;
    }
    else if (i > 0 && arg[0] != '-')
    {
      ExperimentParameters::path_to_level = arg;
    }
  }
}

Parameters Experiment::init_params()
{
  Parameters res;
  
  if (strcmp(ExperimentParameters::param_filename.c_str(), "") != 0) { 
    res.Load(ExperimentParameters::param_filename.c_str());
  }
  
  return res;
}

void Experiment::parse_experiment_parameters()
{
  string line, s;
  double d;
  
  std::ifstream stream(ExperimentParameters::experimentparam_filename);
  
  
  
  while (std::getline(stream, line)) {
    stringstream ss;
    
    ss << line;
    
    if (!(ss >> s >> d)) continue;
    
    if (s == "maxgens") 			ExperimentParameters::max_gens = (int) d;
    else if (s == "numcores")			ExperimentParameters::num_cores = (int) d;
    else if (s == "randseed") {
      ExperimentParameters::using_seed = true;
      ExperimentParameters::seed = (int) d;
    }
    else if (s == "autosaveinterval")		
      if ((int) d == 0)				ExperimentParameters::autosave_best = true;
      else 					ExperimentParameters::autosave_interval = (int) d;
    
    else if (s == "numrangesensors") 		ExperimentParameters::AMOUNT_RANGE_SENSORS = (int) d;
    else if (s == "numdepthsensors") 		ExperimentParameters::AMOUNT_DEPTH_SENSORS = (int) d;
    else if (s == "numpiesensors") 		ExperimentParameters::AMOUNT_PIESLICE_SENSORS = (int) d;
    
    else if (s == "spacingdepthsensors")	ExperimentParameters::SPACING_DEPTH_SENSORS = (int) d;
    else if (s == "radiuspiesensors")		ExperimentParameters::RADIUS_PIESLICE_SENSORS = (int) d;
    
    else if (s == "numhiddenstartneurons")	ExperimentParameters::num_hidden_start_neurons = (int) d;
    
    else if (s == "numhiddenstartneuronscppn")	ExperimentParameters::num_hidden_start_neurons_cppn = (int) d;
    
    else if (s == "noveltysearch" && (int) d)	ExperimentParameters::novelty_search = true;
    else if (s == "hyperneat" && (int) d)	ExperimentParameters::hyperneat = true;
  }
}

void Experiment::save_pop(string filename)
{
  mkdir("popfiles", S_IRWXU);
  
  std::stringstream ss = std::stringstream();
  ss << "popfiles/";
  ss << filename;
    pop.Save(ss.str().c_str());
}

void Experiment::fix_filenames()
{
  try {
    ExperimentParameters::path_to_level = boost::filesystem::canonical(ExperimentParameters::path_to_level).string();
  } catch (const std::exception& e) {
    std::cerr << "Bad level file specified" << std::endl;
    exit(-1);
  }
  
  if (strcmp(ExperimentParameters::pop_filename.c_str(), "")) {
    try {
      ExperimentParameters::pop_filename = boost::filesystem::canonical(ExperimentParameters::pop_filename).string();
    } catch (const std::exception& e) {
      std::cerr << "Bad population file specified" << std::endl;
      exit(-1);
    }
  }
  
  if (strcmp(ExperimentParameters::param_filename.c_str(), "")) {
    try {
      ExperimentParameters::param_filename = boost::filesystem::canonical(ExperimentParameters::param_filename).string();
    } catch (const std::exception& e) {
      std::cerr << "Bad parameter file specified" << std::endl;
      exit(-1);
    }
  }
  
  if (strcmp(ExperimentParameters::experimentparam_filename.c_str(), "")) {
    try {
    ExperimentParameters::experimentparam_filename = boost::filesystem::canonical(ExperimentParameters::experimentparam_filename).string();
    } catch (const std::exception& e) {
      std::cerr << "Bad experiment parameter file specified" << std::endl;
      exit(-1);
    }
  }
}

int Experiment::busy_handler(void *data, int retry)
{
  std::cout << "Busy handler";
  if (retry < max_db_retry) {
    std::cout << ", try #" << retry << std::endl;
    sqlite3_sleep(db_sleeptime);
    return 1;
  } else {
    return 0;
  }
}

void Experiment::print_usage()
{
  stringstream ss;
  ss << std::endl;
  ss << "====================================================" << std::endl;
  ss << "== Help for SuperTux + NEAT experiment executable ==" << std::endl;
  ss << "== Available arguments:                           ==" << std::endl;
  ss << "====================================================" << std::endl;
  ss << "--paramfile\t\tMultiNEAT parameter file" << std::endl;
  ss << "--experimentparamfile\tCustom experiment parameter file. Check README for usage" << std::endl;
  ss << std::endl;
  ss << "The following arguments can be overriden by the experiment parameter file:" << std::endl;
  ss << "--seed\t\t\tSeed for MultiNEAT" << std::endl;
  ss << "--noveltysearch\t\tUse Novelty Search instead of objective-based fitness for evolution" << std::endl;
  
  std::cout << ss.str();
}