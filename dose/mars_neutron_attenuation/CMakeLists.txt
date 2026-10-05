# CMakeLists.txt -- MarsNeutronAttenuation
#
#   mkdir build && cd build
#   cmake -DGeant4_DIR=$HOME/repos/geant4-v11.3.2-install/lib/cmake/Geant4 ..
#   make -j$(sysctl -n hw.ncpu)
#
# (Geant4_DIR must point at the directory containing Geant4Config.cmake --
# lib/cmake/Geant4, NOT lib/Geant4-11.3.2, which is a red herring symlink --
# see the mars_gcr_dose_sim build troubleshooting from earlier this session.)

cmake_minimum_required(VERSION 3.16...3.27)
project(MarsNeutronAttenuation)

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(MarsNeutronAttenuation main.cc ${sources} ${headers})
target_link_libraries(MarsNeutronAttenuation ${Geant4_LIBRARIES})
