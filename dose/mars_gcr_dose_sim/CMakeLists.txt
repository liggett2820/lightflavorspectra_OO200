# CMakeLists.txt -- MarsGCRDoseSim
#
# Standard Geant4 application CMake setup (same structure as Geant4's B1
# example). Build against the user's own Geant4 v11.3.2 install:
#
#   mkdir build && cd build
#   cmake -DGeant4_DIR=$HOME/repos/geant4-v11.3.2-install/lib/Geant4-11.3.2 ..
#   make -j$(sysctl -n hw.ncpu)
#
# (adjust the Geant4-11.3.2 subdirectory name under lib/ if it differs --
# check `ls ~/repos/geant4-v11.3.2-install/lib/` if cmake can't find it.)

cmake_minimum_required(VERSION 3.16...3.27)
project(MarsGCRDoseSim)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(MarsGCRDoseSim main.cc ${sources} ${headers})
target_link_libraries(MarsGCRDoseSim ${Geant4_LIBRARIES})
