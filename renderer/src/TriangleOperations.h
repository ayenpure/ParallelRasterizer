#ifndef TRIANGLEOPERATIONS_H_
#define TRIANGLEOPERATIONS_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Utilities.h"
#include "Camera.h"
#include "MatrixOperations.h"
#include <iostream>

using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

class Triangle {
public:
void Print(std::ostream &o)
{
   o << "Triangle: (" << X[0] << ", " << Y[0] << ", " << Z[0]
                  << "/" << X[1] << ", " << Y[1] << ", " << Z[1]
                  << "/" << X[2] << ", " << Y[2] << ", " << Z[2] << endl;
};
	double X[3];
	double Y[3];
	double Z[3];
	double colors[3][3];
	double normals[3][3];
	double shading[3];
	/*
	 * The offset_vertex is the vertex[X,Y] that does not belong to the base
	 * The left_vertex is the vertex[X,Y] that is to the left side of the base
	 * The right_vertex is the vertex[X,Y] that is to the rigth side of the base
	 * slope_on right will hold the slope of the triangle's side on the right and
	 * similarly slope_on_left will hold the slope of the side on the left of the
	 * triangle
	 * left_offset and right_offset are the constants in the line equation for the
	 * sides of the triangle
	 */
	double offset_vertex[2];
	double left_vertex[2], slope_on_left, left_offset;
	double right_vertex[2], slope_on_right, right_offset;
	int offset_index, left_index, right_index;

	double getlowestY() {
		double min = (Y[1] < Y[2]) ? Y[1] : Y[2];
		return (Y[0] < min) ? Y[0] : min;
	}

	double gethighestY() {
		double max = (Y[1] > Y[2]) ? Y[1] : Y[2];
		return (Y[0] > max) ? Y[0] : max;
	}

	double getlowestX() {
		double min = (X[1] < X[2]) ? X[1] : X[2];
		return (X[0] < min) ? X[0] : min;
	}

	double gethighestX() {
		double max = (X[1] > X[2]) ? X[1] : X[2];
		return (X[0] > max) ? X[0] : max;
	}

	double getlowestZ() {
		double min = (Z[1] < Z[2]) ? Z[1] : Z[2];
		return (Z[0] < min) ? Z[0] : min;
	}

	double gethighestZ() {
		double max = (Z[1] > Z[2]) ? Z[1] : Z[2];
		return (Z[0] > max) ? Z[0] : max;
	}



	void determine_triangle_orientation() {
		/*
		 * We determine the vertices that represent the base and the other
		 * offset vertex here.
		 */
		offset_index = ((Y[0] == Y[1]) ? 2 : ((Y[1] == Y[2]) ? 0 : 1));
		offset_vertex[0] = X[offset_index];
		offset_vertex[1] = Y[offset_index];
		if (X[(offset_index + 1) % 3] < X[(offset_index + 2) % 3]) {
			left_index = (offset_index + 1) % 3;
			right_index = (offset_index + 2) % 3;
		} else {
			left_index = (offset_index + 2) % 3;
			right_index = (offset_index + 1) % 3;
		}
		left_vertex[0] = X[left_index];
		left_vertex[1] = Y[left_index];
		right_vertex[0] = X[right_index];
		right_vertex[1] = Y[right_index];

		/*
		 * We solve for the slope and the offset for the sides of the triangle.
		 * the slope is calculated by two point form slope = (Y2 - Y1) / (X2 - X1)
		 * And the offsets are calculated as c = y - slope * x;
		 */
		double diff;
		diff = (offset_vertex[0] - left_vertex[0]);
		slope_on_left =
				(diff == 0) ? 0 : (offset_vertex[1] - left_vertex[1]) / diff;
		left_offset = left_vertex[1] - (double(slope_on_left * left_vertex[0]));

		diff = (offset_vertex[0] - right_vertex[0]);
		slope_on_right =
				(diff == 0) ? 0 : (offset_vertex[1] - right_vertex[1]) / diff;
		right_offset = right_vertex[1]
				- (double(slope_on_right * right_vertex[0]));

	}

	double get_left_x_intercept(int y) {
		/*
		 * We solve for the intercept for the scan line on the left side of the
		 * triangle
		 */
		double intercept_left =
				(slope_on_left == 0) ?
						getlowestX() :
						(double(y - left_offset)) / slope_on_left;
		return intercept_left;
	}

	double get_right_x_intercept(int y) {
		/*
		 * We solve for the intercept for the scan line on the right side of the
		 * triangle
		 */
		double intercept_right =
				(slope_on_right == 0) ?
						gethighestX() :
						(double(y - right_offset)) / slope_on_right;
		return intercept_right;
	}

	bool is_flat_bottom_triangle() {
		/*
		 * Method to determine if the triangle is a flat bottom or a flat
		 * top triangle
		 */
		if (Y[0] == Y[1] || Y[1] == Y[2] || Y[2] == Y[0])
			return true;
		return false;
	}

	void split_triangle(Triangle *t1, Triangle *t2) {
		/*
		 * Method to split an arbitrary traingle
		 */
		double top_vertex[2];
		double bottom_vertex[2];
		double middle_vertex[2];
		int top_index, bottom_index, middle_index;
		top_index = (Y[0] > Y[1] && Y[0] > Y[2]) ? 0 :
					(Y[1] > Y[0] && Y[1] > Y[2]) ? 1 : 2;
		if (Y[(top_index + 1) % 3] > Y[(top_index + 2) % 3]) {
			middle_index = (top_index + 1) % 3;
			bottom_index = (top_index + 2) % 3;
		} else {
			middle_index = (top_index + 2) % 3;
			bottom_index = (top_index + 1) % 3;
		}
		top_vertex[0] = X[top_index];
		top_vertex[1] = Y[top_index];
		bottom_vertex[0] = X[bottom_index];
		bottom_vertex[1] = Y[bottom_index];
		middle_vertex[0] = X[middle_index];
		middle_vertex[1] = Y[middle_index];

		double diff = top_vertex[0] - bottom_vertex[0];
		double slope =
				(diff == 0) ? 0 : (top_vertex[1] - bottom_vertex[1]) / diff;
		double offset = top_vertex[1] - (double(slope * top_vertex[0]));

		double split_vertex[2];
		split_vertex[1] = middle_vertex[1];
		split_vertex[0] =
				(slope == 0) ?
						top_vertex[0] :
						(double(split_vertex[1] - offset)) / slope;

		/*
		 * Calculate for the value of Z buffer and colors for the new vertex
		 */
		double z_split = interpolate(top_vertex[1], bottom_vertex[1],
				Z[top_index], Z[bottom_index], split_vertex[1]);
		double shade_split = interpolate(top_vertex[1], bottom_vertex[1],
				shading[top_index], shading[bottom_index], split_vertex[1]);

		double split_colors[3];
		split_colors[0] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][0], colors[bottom_index][0], split_vertex[1]);
		split_colors[1] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][1], colors[bottom_index][1], split_vertex[1]);
		split_colors[2] = interpolate(top_vertex[1], bottom_vertex[1],
				colors[top_index][2], colors[bottom_index][2], split_vertex[1]);

		double split_normal[3];
		interpolate_vector(top_vertex[1], bottom_vertex[1], normals[top_index],
				normals[bottom_index], split_vertex[1], split_normal);

		t1->X[0] = top_vertex[0];
		t1->Y[0] = top_vertex[1];
		t1->X[1] = middle_vertex[0];
		t1->Y[1] = middle_vertex[1];
		t1->X[2] = split_vertex[0];
		t1->Y[2] = split_vertex[1];
		t1->Z[0] = Z[top_index];
		t1->Z[1] = Z[middle_index];
		t1->Z[2] = z_split;
		t1->shading[0] = shading[top_index];
		t1->shading[1] = shading[middle_index];
		t1->shading[2] = shade_split;
		/*memcpy(t1->normals[0], this->normals[top_index], 3 * sizeof(double));
		 memcpy(t1->normals[1], this->normals[middle_index], 3 * sizeof(double));
		 memcpy(t1->normals[2], split_normal, 3 * sizeof(double));*/

		t1->colors[0][0] = colors[top_index][0];
		t1->colors[0][1] = colors[top_index][1];
		t1->colors[0][2] = colors[top_index][2];
		t1->colors[1][0] = colors[middle_index][0];
		t1->colors[1][1] = colors[middle_index][1];
		t1->colors[1][2] = colors[middle_index][2];
		t1->colors[2][0] = split_colors[0];
		t1->colors[2][1] = split_colors[1];
		t1->colors[2][2] = split_colors[2];

		t2->X[0] = bottom_vertex[0];
		t2->Y[0] = bottom_vertex[1];
		t2->X[1] = middle_vertex[0];
		t2->Y[1] = middle_vertex[1];
		t2->X[2] = split_vertex[0];
		t2->Y[2] = split_vertex[1];
		t2->Z[0] = Z[bottom_index];
		t2->Z[1] = Z[middle_index];
		t2->Z[2] = z_split;
		t2->shading[0] = shading[bottom_index];
		t2->shading[1] = shading[middle_index];
		t2->shading[2] = shade_split;
		/*memcpy(t2->normals[0], this->normals[bottom_index], 3 * sizeof(double));
		 memcpy(t2->normals[1], this->normals[middle_index], 3 * sizeof(double));
		 memcpy(t2->normals[2], split_normal, 3 * sizeof(double));*/

		t2->colors[0][0] = colors[bottom_index][0];
		t2->colors[0][1] = colors[bottom_index][1];
		t2->colors[0][2] = colors[bottom_index][2];
		t2->colors[1][0] = colors[middle_index][0];
		t2->colors[1][1] = colors[middle_index][1];
		t2->colors[1][2] = colors[middle_index][2];
		t2->colors[2][0] = split_colors[0];
		t2->colors[2][1] = split_colors[1];
		t2->colors[2][2] = split_colors[2];
	}

	/*
	 * This method helps to calculate the colors at extremes of the scanline
	 * corresponding to quest_point
	 */
	void calculate_color_for_scanline_extremes(double quest_point,
			double *color_at_left_intercept, double *color_at_right_intercept) {

		color_at_left_intercept[0] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][0], colors[left_index][0],
				quest_point);
		color_at_left_intercept[1] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][1], colors[left_index][1],
				quest_point);
		color_at_left_intercept[2] = interpolate(offset_vertex[1],
				left_vertex[1], colors[offset_index][2], colors[left_index][2],
				quest_point);

		color_at_right_intercept[0] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][0],
				colors[right_index][0], quest_point);
		color_at_right_intercept[1] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][1],
				colors[right_index][1], quest_point);
		color_at_right_intercept[2] = interpolate(offset_vertex[1],
				right_vertex[1], colors[offset_index][2],
				colors[right_index][2], quest_point);
	}

  /*bool calculate_normal(double *normal) {
		  int i = 0;
			int adj_1 = (i + 1) % 3;
			int adj_2 = (i + 2) % 3;
			double adj_1_vector[3] = { X[adj_1] - X[i], Y[adj_1] - Y[i],
					Z[adj_1] - Z[i] };
			normalize_vector(adj_1_vector);
			double adj_2_vector[3] = { X[adj_2] - X[i], Y[adj_2] - Y[i],
					Z[adj_2] - Z[i] };
			normalize_vector(adj_2_vector);
			cross_product(adj_1_vector, adj_2_vector, normal);
			normalize_vector(normal);
      //std::cout << "{" << normal[0] << ", " << normal[1] << ", " << normal[2] << "}" << std::endl;
      if(std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2]))
        return false;
      return true;
	}*/

	void calculate_normals() {
		//cout << "Caluclating normal ----> calculate_normals" << endl;
		for (int i = 0; i < 3; i++) {
			int adj_1 = (i + 1) % 3;
			int adj_2 = (i + 2) % 3;
			double adj_1_vector[3] = { X[adj_1] - X[i], Y[adj_1] - Y[i],
					Z[adj_1] - Z[i] };
			normalize_vector(adj_1_vector);
			double adj_2_vector[3] = { X[adj_2] - X[i], Y[adj_2] - Y[i],
					Z[adj_2] - Z[i] };
			normalize_vector(adj_2_vector);
			double normal[3];
			cross_product(adj_1_vector, adj_2_vector, normal);
			double t_area = vector_magnitude(normal);
			//cout << "Area of Triangle : " << t_area << endl;
			normalize_vector(normal);
			normals[i][0] = normal[0];
			normals[i][1] = normal[1];
			normals[i][2] = normal[2];
		}
	}

};

std::vector<Triangle> SplitTriangle(std::vector<Triangle> list, int parts);

Triangle rotate_triangle(Triangle t, double angle, char axis);

std::vector<Triangle> get_all_octants(Triangle t);

std::vector<Triangle> GetTriangles(const char *filename, char *variable);

std::vector<Triangle> GetTrianglesFromFiles(int no_of_procs, char *variable);

std::vector<Triangle> GetTriangles(int split_parts, int split_rec, int tri_grain);

std::vector< std::vector<Triangle> > GetTrianglesForProcs(int split_parts, int split_rec, int tri_grain);

std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list);

void transformTriangle(Triangle *t, Matrix composite, Camera camera);

#endif
