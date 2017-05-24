#ifndef CAMERAPOSITIONS_H_
#define CAMERAPOSITIONS_H_

double** get_camera_positions(int config_id, int *num_cameras);

void get_camera_position_and_focus(int config_id, int cam_index,
   int num_cameras ,double** camera_positions,double* camera_position,
   double* focus_point);

#endif
