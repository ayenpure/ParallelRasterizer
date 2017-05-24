#include <iostream>
#include <cmath>
#include <algorithm>
#include <string.h>
#include <vector>
#include <sstream>
#include "Utilities.h"

using std::cout;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

double ceil441(double f) {
	return ceil(f - 0.00001);
}

double floor441(double f) {
	return floor(f + 0.00001);
}

void vector_copy(double *source, double *destination, int num) {
	for(int i=0;i<num;i++) {
		destination[i] = source[i];
	}
}

void rotate(double* camera_position, double angle, char axis,double* rotated_camera) {
  double rotation_matrix[3][3];
  if(axis == 'x') {
    double x_rotation_matrix[3][3] = {
      1,0,0,
      0,cos(angle),-sin(angle),
      0,sin(angle),cos(angle)
    };
    memcpy(rotation_matrix, x_rotation_matrix, 9*sizeof(double));
  } else if (axis == 'y') {
    double y_rotation_matrix[3][3] = {
      cos(angle),0,-sin(angle),
      0,1,0,
      sin(angle),0,cos(angle)
    };
    memcpy(rotation_matrix, y_rotation_matrix, 9*sizeof(double));
  } else if (axis == 'z') {
    double z_rotation_matrix[3][3] = {
      cos(angle),-sin(angle),0,
      sin(angle),cos(angle),0,
      0,0,1
    };
    memcpy(rotation_matrix, z_rotation_matrix, 9*sizeof(double));
  }
  /*for (int i = 0 ; i < 4 ; i++)
  {
      char str[256];
      sprintf(str, "(%.7f %.7f %.7f)\n", rotation_matrix[i][0], rotation_matrix[i][1], rotation_matrix[i][2]);
      cout << str;
  }*/
  for(int i = 0; i < 3; i++) {
    rotated_camera[i] = 0;
    for(int j = 0; j < 3; j++) {
      rotated_camera[i] = rotated_camera[i] + camera_position[j]*rotation_matrix[j][i];
    }
    if(abs(rotated_camera[i]) < 0.0000001)
      rotated_camera[i] = 0;
  }
  //cout << "rotated position {" << rotated_camera[0] << ", " << rotated_camera[1] << ", " << rotated_camera[2] << " }" << endl;
}

/**
 * This function interpolates the value at quest_point based on
 * values value_1 and value_2 at points point_1 and point_2 respectively
 */
double interpolate(double point_1, double point_2, double value_1,
		double value_2, double quest_point) {
	double diff = point_2 - point_1;
	double change_prop =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double value_at_quest = value_1 + (change_prop * (value_2 - value_1));
	return value_at_quest;
}

double vector_magnitude(double* quest_vec) {
	return sqrt( (quest_vec[0] * quest_vec[0])
						 + (quest_vec[1] * quest_vec[1])
						 + (quest_vec[2] * quest_vec[2]));
}

void normalize_vector(double* quest_normal) {
	double norm = sqrt(
			(quest_normal[0] * quest_normal[0])
					+ (quest_normal[1] * quest_normal[1])
					+ (quest_normal[2] * quest_normal[2]));
	quest_normal[0] = quest_normal[0] / norm;
	quest_normal[1] = quest_normal[1] / norm;
	quest_normal[2] = quest_normal[2] / norm;
}

double dot_product(double* vector_1, double* vector_2) {
	return (vector_1[0] * vector_2[0]) + (vector_1[1] * vector_2[1])
			+ (vector_1[2] * vector_2[2]);
}

void cross_product(double* vector_1, double* vector_2, double *cross_vec) {
	double product[3] = { (vector_1[1] * vector_2[2])
			- (vector_1[2] * vector_2[1]), (vector_1[2] * vector_2[0])
			- (vector_1[0] * vector_2[2]), (vector_1[0] * vector_2[1])
			- (vector_1[1] * vector_2[0]) };
	vector_copy(product, cross_vec, 3);
}

void print_vector(double *print_vector) {
	cout << "{" << print_vector[0] << ", " << print_vector[1] << ", "
			<< print_vector[2] << "}" << endl;
}

double cot(double angle) {
	return (1 / tan(angle));
}

/**
 * This function interpolates over the vectors normal_1 and normal_2 to yield
 * the vector quest_normal for the point (quest_point[0],quest_point[1]) in space.
 */
void interpolate_vector(double point_1, double point_2, double* normal_1,
		double* normal_2, double quest_point, double* quest_normal) {
	double diff = point_2 - point_1;
	double proportion =
			(diff != 0) ? (quest_point - point_1) / (point_2 - point_1) : 0;
	double diff_vector[3] = { normal_2[0] - normal_1[0], normal_2[1]
			- normal_1[1], normal_2[2] - normal_1[2] };
	diff_vector[0] = proportion * diff_vector[0];
	diff_vector[1] = proportion * diff_vector[1];
	diff_vector[2] = proportion * diff_vector[2];
	quest_normal[0] = normal_1[0] + diff_vector[0];
	quest_normal[1] = normal_1[1] + diff_vector[1];
	quest_normal[2] = normal_1[2] + diff_vector[2];
	normalize_vector(quest_normal);
}

void getRGBAforDepth(double depth, double *toReturn) {
  double enc[] = {1., 255., 65025., 16581375.};
  double intpart;
  for(int i=0;i<4;i++) {
    enc[i] *= depth;
    enc[i] = modf(enc[i], &intpart);
  }
  vector_copy(enc, toReturn, 4);
}

double range_transform(double r2min, double r2max, double r1min, double r1max, double val) {
	double scale = (r2max - r2min) / (r1max - r1min);
	return (val - r1min)*scale + r2min;
}
