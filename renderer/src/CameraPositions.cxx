#include <iostream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "CameraPositions.h"
#include "Utilities.h"
#include "TriangleOperations.h"

#define MAX_GRAINS 12
#define MAX_LEVELS 3
#define SEED_X_VALUE 10
#define LENGTH 20
#define HEIGHT 20
#define WIDTH 20

using std::sin;
using std::cos;
using std::endl;
using std::cout;
using std::abs;

void get_config_random(double *camera_positions[3]) {
  double calc_camera_positions[2*MAX_LEVELS-1][MAX_GRAINS][3];
  double about_y = 2*M_PI/MAX_GRAINS;
  double about_z = M_PI/(2*MAX_LEVELS);
  double seed_position[3] = {SEED_X_VALUE,0,0};
  int middle_index = MAX_LEVELS - 1;
  vector_copy(seed_position,calc_camera_positions[middle_index][0],3);
  for(int i = 0; i < middle_index; i++) {
    rotate(seed_position, (middle_index-i)*about_z, 'z', calc_camera_positions[i][0]);
  }
  for(int i = middle_index+1; i < 2*MAX_LEVELS-1; i++) {
    rotate(seed_position, 2*M_PI-((i-middle_index)*about_z), 'z', calc_camera_positions[i][0]);
  }
  for(int i=0;i<2*MAX_LEVELS-1;i++) {
    for(int j=1;j<MAX_GRAINS;j++){
      rotate(calc_camera_positions[i][0], j*about_y, 'y', calc_camera_positions[i][j]);
    }
  }
  int total_num_cams = (2*MAX_LEVELS-1)*MAX_GRAINS + 2;
  rotate(seed_position, M_PI/2, 'z', camera_positions[0]);
  rotate(seed_position, 3*M_PI/2, 'z', camera_positions[total_num_cams-1]);
  int position = 1;
  for(int i=0;i<2*MAX_LEVELS-1;i++) {
    for(int j=0;j<MAX_GRAINS;j++){
        vector_copy(calc_camera_positions[i][j],camera_positions[position++],3);
    }
  }
  cout << "";
}

void get_uniform_cinema_distribution(double *camera_positions[3]) {
  std::vector< std::vector<Triangle> > proc_triangles;
  proc_triangles = GetTrianglesForProcs(2,1,0);
  cout << proc_triangles.size() << endl;
  for(int i=0;i<proc_triangles.size();i++) {
    Triangle t = proc_triangles[i][0];
    /*for(int k = 0; k < 3; k++) {
      double vertex[3] = {
        t.X[k], t.Y[k], t.Z[k]
      };
      normalize_vector(vertex);
      t.X[k] = vertex[0];
      t.Y[k] = vertex[1];
      t.Z[k] = vertex[3];
    }*/
    camera_positions[i][0] = (t.X[0] + t.X[1] + t.X[2]) / 3;
    camera_positions[i][1] = (t.Y[0] + t.Y[1] + t.Y[2]) / 3;
    camera_positions[i][2] = (t.Z[0] + t.Z[1] + t.Z[2]) / 3;
    normalize_vector(camera_positions[i]);
    camera_positions[i][0] = camera_positions[i][0] * 40;
    camera_positions[i][1] = camera_positions[i][1] * 40;
    camera_positions[i][2] = camera_positions[i][2] * 40;
  }
}

void get_config_helix(double *camera_positions[3]) {
  double y = -15;
  int cam_index = 0;
  double camera_position[3] = {40, y,0};
  double rotated_camera[3];
  camera_positions[cam_index][0] = camera_position[0];
  camera_positions[cam_index][1] = camera_position[1];
  camera_positions[cam_index][2] = camera_position[2];
  //cout << camera_position[0] << ", " << camera_position[1] << ", " << camera_position[2] << endl;
  cam_index++;
  do {
    rotate(camera_position, M_PI / 8, 'y', rotated_camera);
    camera_position[0] = rotated_camera[0];
    y += .375;
    camera_position[1] = y;
    camera_position[2] = rotated_camera[2];
    //cout << camera_position[0] << ", " << camera_position[1] << ", " << camera_position[2] << endl;
    camera_positions[cam_index][0] = camera_position[0];
    camera_positions[cam_index][1] = camera_position[1];
    camera_positions[cam_index][2] = camera_position[2];
    cam_index++;
  } while ( y < 15);
}

void get_camera_positions_for_shawn(double *camera_positions[3]) {
  camera_positions[0][0] = 0;
  camera_positions[0][1] = 0;
  camera_positions[0][2] = WIDTH / 2.;
  camera_positions[1][0] = 0;
  camera_positions[1][1] = 0;
  camera_positions[1][2] = - WIDTH / 2.;
  camera_positions[2][0] = 0;
  camera_positions[2][1] = HEIGHT / 2.;
  camera_positions[2][2] = 0;
  camera_positions[3][0] = 0;
  camera_positions[3][1] = - HEIGHT / 2.;
  camera_positions[3][2] = 0;
  camera_positions[4][0] = LENGTH / 2.;
  camera_positions[4][1] = 0;
  camera_positions[4][2] = 0;
  camera_positions[5][0] = - LENGTH / 2.;
  camera_positions[5][1] = 0;
  camera_positions[5][2] = 0;
}

void get_camera_position_and_focus(int config_id, int cam_index, int num_cameras ,double** camera_positions,double* camera_position,double* focus_point) {
  int focus_index = 0;
  camera_position[0] = camera_positions[cam_index][0];
  camera_position[1] = camera_positions[cam_index][1];
  camera_position[2] = camera_positions[cam_index][2];
  switch(config_id) {
    default :
    case 0 :
      focus_point[0] = 0;
      focus_point[1] = 0;
      focus_point[2] = 0;
    break;
    case 1 :
      if(cam_index < 16)
        focus_index = (cam_index + 6) % 16;
      else {
        int quotient = (cam_index - 16) / 14;
        int pseudo_index = (cam_index - 16) % 14;
        int pseudo_focus_index = (pseudo_index + 5) % 14;
        focus_index = ((quotient*14) + 16) + pseudo_focus_index;
      }
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 2 :
      if(cam_index < 16)
        focus_index = (cam_index + 7) % 16;
      else {
        int quotient = (cam_index - 16) / 14;
        int pseudo_index = (cam_index - 16) % 14;
        int pseudo_focus_index = (pseudo_index + 6) % 14;
        focus_index = ((quotient*14) + 16) + pseudo_focus_index;
      }
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 3 :
      focus_index = cam_index + 7;
      focus_point[0] = camera_positions[focus_index][0];
      focus_point[1] = camera_positions[focus_index][1];
      focus_point[2] = camera_positions[focus_index][2];
    break;
    case 4 :
      focus_index = cam_index + 7;
      focus_point[0] = 0;
      focus_point[1] = camera_positions[cam_index][1];
      focus_point[2] = 0;
    break;
    case 6:
    focus_point[0] = camera_positions[cam_index][0];
    focus_point[1] = camera_positions[cam_index][1];
    focus_point[2] = camera_positions[cam_index][2];
    camera_position[0] = 0;
    camera_position[1] = 0;
    camera_position[2] = 0;
    break;
  }
}

double** allocate_memory(int num_cameras) {
  double **camera_positions;
  camera_positions = (double**)malloc(num_cameras*sizeof(double*));
  for(int i = 0; i < num_cameras; i++)
    camera_positions[i] = (double*)malloc(3*sizeof(double));
  return camera_positions;
}


double** get_camera_positions(int config_id, int *num_cameras) {
  double **camera_positions;
  switch(config_id) {
    case 0:
    case 1:
    case 2:
    default :
      *num_cameras = (2*MAX_LEVELS-1)*MAX_GRAINS + 2;
      camera_positions = allocate_memory(*num_cameras);
      get_config_random(camera_positions);
    break;
    case 3:
    case 4:
      *num_cameras = 74;
      camera_positions = allocate_memory(81);
      get_config_helix(camera_positions);
    break;
    case 5:
      *num_cameras = 64;
      camera_positions = allocate_memory(64);
      get_uniform_cinema_distribution(camera_positions);
    break;
    case 6:
      *num_cameras = 6;
      camera_positions = allocate_memory(6);
      get_camera_positions_for_shawn(camera_positions);
    break;
  }
  return camera_positions;
}
