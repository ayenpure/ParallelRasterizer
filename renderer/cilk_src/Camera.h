#ifndef CAMERA_H_
#define CAMERA_H_

#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>
#include "MatrixOperations.h"
#include "Screen.h"

class Camera {
public:
	double near, far;
	double angle;
	double position[3];
	double focus[3];
	double up[3];

	Matrix CameraTransform(void) {
		double z_vector[3] = { position[0] - focus[0], position[1] - focus[1],
				position[2] - focus[2] };
		if( (position[0] == 0 && position[2] == 0 ) && (focus[0] == 0 && focus[2] == 0 )) {
			if(position[1] < focus[1]) {
					up[0] = 0; up[1] = 0; up[2] = 1;
			} else {
					up[0] = 0; up[1] = 0; up[2] = -1;
			}
		}
		normalize_vector(z_vector);
		double x_vector[3], y_vector[3];
		cross_product(up, z_vector, x_vector);
		normalize_vector(x_vector);
		cross_product(z_vector, x_vector, y_vector);
		normalize_vector(y_vector);
		double frame_vec[3] = { 0 - position[0], 0 - position[1], 0
				- position[2] };
		//normalize_vector(frame_vec);
		/*cout << "Camera Position (O) :";
		print_vector(position);
		cout << "Focus (F) :";
		print_vector(focus);
		cout << "V1 (U) :";
		print_vector(x_vector);
		cout << "V2 (V) :";
		print_vector(y_vector);
		cout << "V3 (W) :";
		print_vector(z_vector);
		cout << "T :";
		print_vector(frame_vec);*/
		double cartesian_x[3] = { 1, 0, 0 };
		double cartesian_y[3] = { 0, 1, 0 };
		double cartesian_z[3] = { 0, 0, 1 };
		double camera_transform[4][4] = { dot_product(x_vector, cartesian_x),
				dot_product(y_vector, cartesian_x), dot_product(z_vector,
						cartesian_x), 0, dot_product(x_vector, cartesian_y),
				dot_product(y_vector, cartesian_y), dot_product(z_vector,
						cartesian_y), 0, dot_product(x_vector, cartesian_z),
				dot_product(y_vector, cartesian_z), dot_product(z_vector,
						cartesian_z), 0, dot_product(x_vector, frame_vec),
				dot_product(y_vector, frame_vec), dot_product(z_vector,
						frame_vec), 1 };
		Matrix m;
		memcpy(m.A, camera_transform, 16 * sizeof(double));
		return m;
	}

	Matrix ViewTransform(void) {
		double view_transform[4][4] = { cot(angle / 2), 0, 0, 0, 0, cot(
				angle / 2), 0, 0, 0, 0, (far + near) / (far - near), -1, 0, 0,
				(2 * far * near) / (far - near), 0 };
		Matrix m;
		memcpy(m.A, view_transform, 16 * sizeof(double));
		return m;
	}

	Matrix DeviceTransform(Screen screen) {
		double device_transform[4][4] = { screen.width / 2., 0, 0, 0, 0,
				screen.height / 2., 0, 0, 0, 0, 1, 0, screen.width / 2.,
				screen.height / 2., 0, 1 };
		Matrix m;
		memcpy(m.A, device_transform, 16 * sizeof(double));
		return m;
	}
	;
};

Matrix get_total_transform_matrix(Matrix camera_transform,
		Matrix view_transform, Matrix device_transform);

Camera GetCamera(double* camera_position, double* focus_point);

#endif
