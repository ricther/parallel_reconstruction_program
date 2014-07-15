##sudo apt-get install libxmu-dev libxi-dev freeglut

PROJECT(APP)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_BUILD_TYPE Debug)
#The version number.
#set(APP_VERSION_MAJOR 1)
#set(APP_VERSION_MINOR 0)

#to the source code 
#configure_file(
#  "${PROJECT_SOURCE_DIR}/code/main.h.in"
#  "${PROJECT_BINARY_DIR}/code/main.h"
#)


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include(FindOpenMP)
if(OPENMP_FOUND)
MESSAGE(FOUND_OPENMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#add the binary tree to the search path for include files
#so we can find the main.h and **.h
include_directories(${APP_SOURCE_DIR}/code)

FILE(GLOB SOURCE_1 "${APP_SOURCE_DIR}/code/*.cpp")


add_executable(APP ${SOURCE_1})

if(VTK_LIBRARIES)
  target_link_libraries(APP ${VTK_LIBRARIES})
else()
  target_link_libraries(APP vtkHybrid)
endif()

Message("FINISHED")
