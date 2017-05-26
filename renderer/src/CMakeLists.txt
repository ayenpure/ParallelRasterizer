cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1F)
SET(VTK_DIR /path/to/VTK6.0.0)
SET(CMAKE_CXX_FLAGS "--std=c++11  ${CMAKE_CXX_FLAGS}")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(renderer Renderer.cxx TriangleOperations.h TriangleOperations.cxx Utilities.h Utilities.cxx MatrixOperations.h
               Screen.h Camera.h Camera.cxx RenderFunctions.h RenderFunctions.cxx LightingParameters.h)

if(VTK_LIBRARIES)
  target_link_libraries(renderer ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkHybrid)
endif()
