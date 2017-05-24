#include "Camera.h"

#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>
#include "MatrixOperations.h"

double SineParameterize(int curFrame, int nFrames, int ramp) {
	int nNonRamp = nFrames - 2 * ramp;
	double height = 1. / (nNonRamp + 4 * ramp / M_PI);
	if (curFrame < ramp) {
		double factor = 2 * height * ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double) curFrame) / ramp);
		return (1. - eval) * factor;
	} else if (curFrame > nFrames - ramp) {
		int amount_left = nFrames - curFrame;
		double factor = 2 * height * ramp / M_PI;
		double eval = cos(M_PI / 2 * ((double) amount_left / ramp));
		return 1. - (1 - eval) * factor;
	}
	double amount_in_quad = ((double) curFrame - ramp);
	double quad_part = amount_in_quad * height;
	double curve_part = height * (2 * ramp) / M_PI;
	return quad_part + curve_part;
}

Camera GetCamera(double* camera_position, double* focus_point) {
	//double t = SineParameterize(frame, nframes, nframes / 10);
	Camera c;
	c.near = 0.01;
	c.far =  300.;
	c.angle = M_PI / 2;
	c.position[0] = camera_position[0];
	c.position[1] = camera_position[1];
	c.position[2] = camera_position[2];
	c.focus[0] = focus_point[0];
	c.focus[1] = focus_point[1];
	c.focus[2] = focus_point[2];
	c.up[0] = 0;
	c.up[1] = 1;
	c.up[2] = 0;
	return c;
}

Matrix get_total_transform_matrix(Matrix camera_transform,
		Matrix view_transform, Matrix device_transform) {
	Matrix composite;
	composite = Matrix::ComposeMatrices(camera_transform, view_transform);
	composite = Matrix::ComposeMatrices(composite, device_transform);
	return composite;
}
