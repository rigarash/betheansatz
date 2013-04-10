#  Copyright Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp> 2013.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)

# The minimum version is stricter than ALPS itself,
# in order to load the correct FindBoost.cmake when testing.
cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)
project(betheansatz NONE)

######################################################################
# find ALPS Library
######################################################################
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

######################################################################
# applications
######################################################################
# enable C and C++ compilers
enable_language(CXX)

# betheansatz
add_executable(betheansatz main.cc)
target_link_libraries(betheansatz ${ALPS_LIBRARIES})

######################################################################
# unit test (using Boost.Test)
######################################################################

# Check and build Boost.Test Unit Test Framework
# since ALPS does not build it by default.
if (NOT ALPS_HAS_BOOST_TEST)
  find_package(Boost 1.48.0 COMPONENTS unit_test_framework)
endif (NOT ALPS_HAS_BOOST_TEST)

if (ALPS_HAS_BOOST_TEST OR Boost_FOUND)
  add_definitions(-DBOOST_TEST_DYN_LINK)
  # If ALPS library is linked with Boost.test Unit Test Framework
  # or precompiled Boost.test Unit Test Framework is found,
  # set ALPS_HAS_BOOST_TEST definition and
  # use header-only version of Boost.Test Unit Test Framework.
  add_definitions(-DALPS_LINK_BOOST_TEST)
endif (ALPS_HAS_BOOST_TEST OR Boost_FOUND)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})
# test betheansatz.h
add_executable(test_betheansatz test_betheansatz.cc)
target_link_libraries(test_betheansatz ${ALPS_LIBRARIES} ${Boost_LIBRARIES})
enable_testing()
add_test(test_betheansatz test_betheansatz)
