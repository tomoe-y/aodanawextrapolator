# Set up the project.
cmake_minimum_required( VERSION 3.10 )
project( MyAnalysisProject VERSION 1.0.0 LANGUAGES CXX )

# Try to figure out what project is our parent. Just using a hard-coded list
# of possible project names. Basically the names of all the other
# sub-directories inside the Projects/ directory in the repository.
set( _parentProjectNames Athena AthenaP1 AnalysisBase AthAnalysis
   AthSimulation AthDerivation AnalysisTop )
set( _defaultParentProject AnalysisBase )
foreach( _pp ${_parentProjectNames} )
   if( NOT "$ENV{${_pp}_DIR}" STREQUAL "" )
      set( _defaultParentProject ${_pp} )
      break()
   endif()
endforeach()

# Set the parent project name based on the previous findings:
set( ATLAS_PROJECT ${_defaultParentProject}
   CACHE STRING "The name of the parent project to build against" )

# Clean up:
unset( _parentProjectNames )
unset( _defaultParentProject )

# Find the base project that my project uses.
find_package( ${ATLAS_PROJECT} REQUIRED )

# Set up the usage of CTest in the project.
atlas_ctest_setup()

# Set up the "ATLAS project" to build all of the packages in the project's source directory.
atlas_project( MyAnalysisProject 1.0.0
   USE ${ATLAS_PROJECT} ${${ATLAS_PROJECT}_VERSION} )

# Generate an environment setup script for the built project.
lcg_generate_env( SH_FILE "${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh" )
install( FILES "${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh"
   DESTINATION . )

# Set up the usage of CPack in the project.
atlas_cpack_setup()
