# Boost including parts taken from MultiNEAT CMakeLists.txt

cmake_minimum_required(VERSION 2.6)
project(stneatexperiment)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

set(Boost_USE_STATIC_LIBS       OFF) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost COMPONENTS
	date_time
	system
	filesystem
	serialization)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(sqlite3)
if (SQLITE3_FOUND)
	include_directories(${SQLITE3_INCLUDE_DIRS})
endif()

include_directories(${root_SOURCE_DIR})
include_directories(multineat)
include_directories(multineat/sensors)


message( "CMAKE_CXX_FLAGS IS ${CMAKE_CXX_FLAGS}" )

add_executable(stneatexperiment multineat/Substrate.cpp multineat/Parameters.cpp multineat/Utils.cpp multineat/NeuralNetwork.cpp multineat/PhenotypeBehavior.cpp multineat/Genome.cpp multineat/Random.cpp multineat/Innovation.cpp multineat/Species.cpp multineat/Population.cpp experiment.cpp multineat/Behavior.cpp multineat/sensors/depthfindersensor.cpp multineat/sensors/pieslicesensor.cpp multineat/sensors/pieslicesensorspecial.cpp multineat/sensors/rangefindersensor.cpp multineat/sensors/sensor.cpp multineat/sensors/sensormanager.cpp parameters.cpp)

target_link_libraries(stneatexperiment ${Boost_LIBRARIES} ${SQLITE3_LIBRARIES})

install(TARGETS stneatexperiment RUNTIME DESTINATION bin)
