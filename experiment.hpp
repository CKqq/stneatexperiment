#ifndef EXPERIMENT_HEADER
#define EXPERIMENT_HEADER

#include <string>
#include <vector>

#include "multineat/Population.h"
#include "multineat/Genome.h"
#include "multineat/Species.h"
#include "multineat/Parameters.h"
#include "multineat/Behavior.h"

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
  static int select_handler_ns(void* data, int argc, char** argv, char** colNames);
  
  static double sparseness(Genome* g);
  void update_sparsenesses();
  
  static int cur_gen;
  
  static int top_genome_id;
  static double top_fitness;
  
  static bool ns;
  
private:  
  // Static for db handler access
  static vector<Genome*> genomes;
  
  // Archive for Novelty Search
  static vector<Genome*> archive;

  void refresh_genome_list();
  
  void save_pop(string filename);
  
  void start_processes();
  void update_db();
  void set_fitness_values();
  
  Genome start_genome;
  Population pop;
  
  static void print_usage();

public:  
static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
  {
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0) {
	fprintf(out, "  <empty, possibly corrupt>\n");
	return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
	char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

	// find parentheses and +address offset surrounding the mangled name:
	// ./module(function+0x15c) [0x8048a6d]
	for (char *p = symbollist[i]; *p; ++p)
	{
	    if (*p == '(')
		begin_name = p;
	    else if (*p == '+')
		begin_offset = p;
	    else if (*p == ')' && begin_offset) {
		end_offset = p;
		break;
	    }
	}

	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';

	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():

	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
					    funcname, &funcnamesize, &status);
	    if (status == 0) {
		funcname = ret; // use possibly realloc()-ed string
		fprintf(out, "  %s : %s+%s\n",
			symbollist[i], funcname, begin_offset);
	    }
	    else {
		// demangling failed. Output function name as a C function with
		// no arguments.
		fprintf(out, "  %s : %s()+%s\n",
			symbollist[i], begin_name, begin_offset);
	    }
	}
	else
	{
	    // couldn't parse the line? print the whole line.
	    fprintf(out, "  %s\n", symbollist[i]);
	}
    }

    free(funcname);
    free(symbollist);
  }
  
static void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  print_stacktrace();
  exit(1);
}
};

#endif