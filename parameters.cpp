#include "parameters.hpp"

string ExperimentParameters::path_to_level = "./bin/data/levels/world1/01\ -\ Welcome\ to\ Antarctica.stl";


int ExperimentParameters::TILE_WIDTH = 32;

int ExperimentParameters::AMOUNT_RANGE_SENSORS = 5;

int ExperimentParameters::AMOUNT_DEPTH_SENSORS = 7;
int ExperimentParameters::SPACING_DEPTH_SENSORS = 64;

int ExperimentParameters::AMOUNT_PIESLICE_SENSORS = 5;

int ExperimentParameters::RADIUS_PIESLICE_SENSORS = 256;

int ExperimentParameters::RADIUS_PIESLICE_SENSORS_SPECIAL = 320;

string ExperimentParameters::pop_filename = "";
string ExperimentParameters::param_filename = "";
string ExperimentParameters::experimentparam_filename = "";

int ExperimentParameters::num_cores = 3;
int ExperimentParameters::max_gens = 100;

bool ExperimentParameters::using_seed = false;
int ExperimentParameters::seed = 0;

int ExperimentParameters::autosave_interval = 0;
bool ExperimentParameters::autosave_best = false;

int ExperimentParameters::num_range_sensors = 5;
int ExperimentParameters::num_depth_sensors = 7;
int ExperimentParameters::num_pieslice_sensors = 7;

int ExperimentParameters::num_hidden_start_neurons = 0;

int ExperimentParameters::num_hidden_start_neurons_cppn = 0;

bool ExperimentParameters::novelty_search = false;
bool ExperimentParameters::hyperneat = false;