cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1F)
SET(VTK_DIR /path/to/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(qualitywork quality_work.cxx CameraPositions.h CameraPositions.cxx)
add_executable(qualitycamera quality_camera.cxx CameraPositions.h CameraPositions.cxx)
add_executable(qualitytriangles quality_triangle.cxx CameraPositions.h CameraPositions.cxx)
add_executable(qualitympiwork quality_mpi_work.cxx CameraPositions.h CameraPositions.cxx)
add_executable(testrender test_renders.cxx CameraPositions.h CameraPositions.cxx)

#if(MPI_COMPILE_FLAGS)
#  set_target_properties(qualitympiwork PROPERTIES
#    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#endif()
#
#if(MPI_LINK_FLAGS)
#  set_target_properties(qualitympiwork PROPERTIES
#    LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()

if(VTK_LIBRARIES)
  target_link_libraries(qualitywork ${VTK_LIBRARIES})
  target_link_libraries(qualitycamera ${VTK_LIBRARIES})
  target_link_libraries(qualitytriangles ${VTK_LIBRARIES})
  target_link_libraries(testrender ${VTK_LIBRARIES})
  if(MPI_C_LIBRARIES)
    target_link_libraries(qualitympiwork ${VTK_LIBRARIES} ${MPI_C_LIBRARIES})
  endif()
else()
  target_link_libraries( vtkHybrid)
endif()
