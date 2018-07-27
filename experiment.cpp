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

#include <sqlite3.h>

#include "experiment.hpp"

#include <signal.h>

string Experiment::temp_pop_filename = "temp_pop";

string Experiment::db_filename = "neat.db";

int Experiment::max_db_retry = 100;
int Experiment::db_sleeptime = 50;

string Experiment::pop_filename = "";
string Experiment::param_filename = "";
string Experiment::experimentparam_filename = "";

string Experiment::path_to_supertux_executable = "./bin/build/supertux2";
string Experiment::path_to_level = "./bin/data/levels/world1/01\ -\ Welcome\ to\ Antarctica.stl";
int Experiment::num_cores = 3;
int Experiment::max_gens = 100;

bool Experiment::using_seed = false;
int Experiment::seed = 0;

int Experiment::cur_gen = 0;
int Experiment::top_genome_id = 0;
double Experiment::top_fitness = 0;

int Experiment::autosave_interval = 0;

int Experiment::num_range_sensors = 5;
int Experiment::num_depth_sensors = 7;
int Experiment::num_pieslice_sensors = 7;

int Experiment::num_hidden_start_neurons = 0;

bool Experiment::ns = false;

Parameters Experiment::params;

vector<Genome*> Experiment::genomes;
vector<PhenotypeBehavior> Experiment::archive;

int main(int argc, char **argv) {  
  signal(SIGSEGV, Experiment::handler);
  
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
start_genome(0, num_range_sensors + num_depth_sensors + num_pieslice_sensors + 1, num_hidden_start_neurons, 
  6, false, UNSIGNED_SIGMOID, UNSIGNED_SIGMOID, 1, params),
  pop(strcmp(Experiment::pop_filename.c_str(), "") ? 
  Population(Experiment::pop_filename.c_str()) : Population(start_genome, params, true, 2.0, (using_seed ? seed : (int) time(0))))
{
}

Experiment::~Experiment()
{
}

void Experiment::run_evolution() {
  std::remove(db_filename.c_str());
  
  for (int i = 1; i <= max_gens; i++) {    
    std::cout << "Working on generation #" << i << "..." << std::endl;
    
    cur_gen = i;
    refresh_genome_list();
    
    if (ns) {
      std::vector<PhenotypeBehavior>* behaviors = new std::vector<PhenotypeBehavior>();
      
      for (int i = 0; i < genomes.size(); i++) 
	behaviors->push_back(Behavior());
      
      pop.InitPhenotypeBehaviorData(behaviors, &archive);
    }
    
    if (autosave_interval != 0 && i % autosave_interval == 0) {
      std::ostringstream ss;
      ss << "./neat_gen" << i;
      save_pop(ss.str().c_str());
    }
    
    // Prepare db file...
    update_db();
    
    // Run supertux processes...
    start_processes();
    
    // ... and get the fitness values from the db file
    set_fitness_values();
    
    std::cout << "Done. Top fitness: " << top_fitness << " by individual #" << top_genome_id << std::endl;
    
    pop.Epoch();
  }
}

void Experiment::start_processes()
{
  // Save pop to single file for all child processes to read

  pop.Save(temp_pop_filename.c_str());
  
  // Distribute genomes as evenly as possible
  int remaining_genome_count = genomes.size();
  int cores_left = num_cores;
  
  std::cout << "Genomes: " << remaining_genome_count << std::endl;
  
  std::vector<int> startindices;
  std::vector<int> endindices;
  
  int curindex = 0;
  
  for (int j = 0; j < num_cores; j++) {
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
  
  // All child processes read from the same population file
  ss << " --popfile " << boost::filesystem::canonical(temp_pop_filename);
  
  // Only add experiment and MultiNEAT param files if they're configured...
  if (strcmp(experimentparam_filename.c_str(), "")) 	ss << " --experimentparamfile " << experimentparam_filename;
  if (strcmp(param_filename.c_str(), "")) 		ss << " --paramfile " << param_filename;
  
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
  
  ss << path_to_level;
  
  ss << "'";
  
//   std::cout << ss.str() << std::endl;
  
  // Ready to rumble
  std::system(ss.str().c_str());

  // Remove temporary population file. If we have to, we'll make a new one
  std::remove(temp_pop_filename.c_str());
}

void Experiment::update_db()
{
  sqlite3* db;  
  sqlite3_open("neat.db", &db);
  
  sqlite3_busy_handler(db, busy_handler, (void*) nullptr);
    
  char* err;
  std::stringstream ss;
  ss << "CREATE TABLE GEN" << cur_gen << "(id INT PRIMARY KEY NOT NULL, fitness REAL, qLeft REAL, qRight REAL, qUp REAL, qDown REAL, qJump REAL, qAction REAL);";
    
  sqlite3_exec(db, ss.str().c_str(), 0, 0, &err);
  
  for(std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    ss.str("");
    
    // Negative fitness if not evaluated yet. If fitness values stay negative, NEAT might just break
    // Last values are for ns behavior
    ss << "INSERT INTO GEN" << cur_gen << " VALUES(" << (*it)->GetID() << ", -1, -1, -1, -1, -1, -1, -1);";
    
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
  
  sqlite3_exec(db, ss.str().c_str(), (ns ? select_handler_ns : select_handler), 0, &err);
  
  if (ns) {
    for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
      if ((*it)->m_PhenotypeBehavior->m_Data[0].size() == 0) {
	std::cout << "ERROR OCCURED WITH GENOME #" << (*it)->GetID() << std::endl;
      }
    }
    update_sparsenesses();
  }

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
//     std::cout << "Updating sparseness of genome #" << (*it)->GetID() << "..." << std::endl;
//     double sparseness = Experiment::sparseness(*it);
    
    (*it)->SetFitness(pop.ComputeSparseness(*(*it)));
    (*it)->SetEvaluated();
  }
}

int Experiment::select_handler(void* data, int argc, char** argv, char** colNames)
{
  int id = std::stoi(argv[0]);
  double fitness = std::stod(argv[1]);
  
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    if ((*it)->GetID() == id) {
      (*it)->SetFitness((Experiment::ns) ? sparseness(*it) : fitness);
      (*it)->SetEvaluated();
      if (fitness > top_fitness) {
	top_fitness = fitness;
	top_genome_id = id;
      }
      
      break;
    }
  }
  
  return 0;
}

int Experiment::select_handler_ns(void* data, int argc, char** argv, char** colNames)
{
  int id = std::stoi(argv[0]);

  double fitness = std::stod(argv[1]);
  
  for (std::vector<Genome*>::iterator it = genomes.begin(); it != genomes.end(); ++it) {
    if ((*it)->GetID() == id) {
            
      for (int i = 2; i < 8; i++) {
	(*it)->m_PhenotypeBehavior->m_Data[0][std::stod(argv[i])];
      }
	
      if (fitness > top_fitness) {
	top_fitness = fitness;
	top_genome_id = id;
      }
      
      break;
    }
  }
  
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
	experimentparam_filename = argv[++i];
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
	param_filename = argv[++i];
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
	pop_filename = argv[++i];
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
	Experiment::using_seed = true;
      } else {
	std::cout << "Please specify a proper seed after '--seed'" << std::endl;
      }
    }
    else if (arg == "--noveltysearch")
    {
      Experiment::ns = true;
    }
  }
}

Parameters Experiment::init_params()
{
  Parameters res;
  
  if (strcmp(param_filename.c_str(), "") != 0) { 
    res.Load(param_filename.c_str());
  }
  
  return res;
}

void Experiment::parse_experiment_parameters()
{
  string line, s;
  double d;
  
  std::ifstream stream(experimentparam_filename);
  
  
  
  while (std::getline(stream, line)) {
    stringstream ss;
    
    ss << line;
    
    if (!(ss >> s >> d)) continue;
    
    if (s == "maxgens") 			Experiment::max_gens = (int) d;
    else if (s == "numcores")			Experiment::num_cores = (int) d;
    else if (s == "randseed") {
      using_seed = true;
      seed = (int) d;
    }
    else if (s == "autosaveinterval")		Experiment::autosave_interval = (int) d;
    
    else if (s == "numrangesensors") 		Experiment::num_range_sensors = (int) d;
    else if (s == "numdepthsensors") 		Experiment::num_depth_sensors = (int) d;
    else if (s == "numpiesensors") 		Experiment::num_pieslice_sensors = (int) d;
    
    else if (s == "numhiddenstartneurons")	Experiment::num_hidden_start_neurons = (int) d;
    
    else if (s == "noveltysearch" && (int) d)	Experiment::ns = true;
  }
}

void Experiment::save_pop(string filename)
{
    pop.Save(filename.c_str());
}

void Experiment::fix_filenames()
{
  try {
    path_to_level = boost::filesystem::canonical(path_to_level).string();
  } catch (const std::exception& e) {
    std::cerr << "Bad level file specified" << std::endl;
    exit(-1);
  }
  
  if (strcmp(pop_filename.c_str(), "")) {
    try {
      pop_filename = boost::filesystem::canonical(pop_filename).string();
    } catch (const std::exception& e) {
      std::cerr << "Bad population file specified" << std::endl;
      exit(-1);
    }
  }
  
  if (strcmp(param_filename.c_str(), "")) {
    try {
      param_filename = boost::filesystem::canonical(param_filename).string();
    } catch (const std::exception& e) {
      std::cerr << "Bad parameter file specified" << std::endl;
      exit(-1);
    }
  }
  
  if (strcmp(experimentparam_filename.c_str(), "")) {
    try {
    experimentparam_filename = boost::filesystem::canonical(experimentparam_filename).string();
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


