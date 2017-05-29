#include <cmath>
#include "Utilities.h"
#include "RenderFunctions.h"
#include "LightingParameters.h"
#include "TriangleOperations.h"

int numFrags = 0;
double calculate_for_specular_lighting(LightingParameters lp, double *V,
		double* N) {
	/* R = 2*(L . N)*N - L
	 * return ==> V . R
	 */
	double two_L_dot_N = 2 * dot_product(lp.lightDir, N);
	double two_L_dot_N_N[3] = { two_L_dot_N * N[0], two_L_dot_N * N[1],
			two_L_dot_N * N[2] };
	double R[3] = { two_L_dot_N_N[0] - lp.lightDir[0], two_L_dot_N_N[1]
			- lp.lightDir[1], two_L_dot_N_N[2] - lp.lightDir[2] };
	return dot_product(R, V);
}

double calculate_phong_shading(LightingParameters lp, double *view_direction,
		double *normal) {
	double diffuse_component = dot_product(lp.lightDir, normal);
	double specular_component = pow(
			calculate_for_specular_lighting(lp, view_direction, normal),
			lp.alpha);
	if (specular_component < 0 || specular_component != specular_component) {
		specular_component = 0;
	}
	double shading_amount = lp.Ka + lp.Kd * abs(diffuse_component)
			+ lp.Ks * specular_component;
}

void scan_line(Triangle *t, Screen *s) {
	/*
	 * Execute scan line algorithm for the current triangle.
	 */
	double pixels_colored = 0;
	double y_min = t->getlowestY();
	double y_max = t->gethighestY();
	double x_min = t->getlowestX();
	double x_max = t->gethighestX();
	double z_min = t->getlowestZ();
	double z_max = t->gethighestZ();

	if (y_max < 0 || y_min > s->height)
		return;
	if (x_max < 0 || x_min > s->width)
		return;
	if (!((z_min > -1 && z_min < 1) && (z_max > -1 && z_max < 1)))
		return;

	if (y_min < 0)
		y_min = 0;
	if (y_max > s->height)
		y_max = s->height;

	if (x_min < 0)
		x_min = 0;
	if (x_max > s->width)
		x_max = s->width;

        double left_intercept;
        double right_intercept;
        double z_left_intercept;
        double z_right_intercept;
        double color_at_left_intercept[3], color_at_right_intercept[3];
        double shading_left_intercept ;
        double shading_right_intercept;
        int current_x;
        double current_z;
        double current_shading;
        // = { 0, 0, 0 };

	// Determine the orientation for the triangle
	t->determine_triangle_orientation();
	// Color the pixels that are inside the triangle
	for (int current_y = ceil441(y_min); current_y <= floor441(y_max);
			current_y++) {

		left_intercept    = t->get_left_x_intercept(current_y);
		right_intercept   = t->get_right_x_intercept(current_y);
		z_left_intercept  = interpolate(t->offset_vertex[1],
				    t->left_vertex[1], t->Z[t->offset_index],
                                    t->Z[t->left_index],
				    current_y);
		z_right_intercept = interpolate(t->offset_vertex[1],
				t->right_vertex[1], t->Z[t->offset_index],
                                t->Z[t->right_index],
				current_y);

		t->calculate_color_for_scanline_extremes(current_y,
				color_at_left_intercept,
                                color_at_right_intercept);

		shading_left_intercept = interpolate(t->offset_vertex[1],
				t->left_vertex[1], t->shading[t->offset_index],
				t->shading[t->left_index], current_y);
		shading_right_intercept = interpolate(t->offset_vertex[1],
				t->right_vertex[1], t->shading[t->offset_index],
				t->shading[t->right_index], current_y);

		for (current_x = ceil441(left_intercept);
				current_x <= floor441(right_intercept);
                                current_x++) {
			 current_z       = interpolate(left_intercept,
                                        right_intercept,
					z_left_intercept, z_right_intercept,
                                        current_x);
			 current_shading = interpolate(left_intercept,
					right_intercept, shading_left_intercept,
					shading_right_intercept, current_x);

			double color_for_current_pixel[3];
			s->calculate_color_for_pixel(left_intercept, right_intercept,
					current_x, color_at_left_intercept,
					color_at_right_intercept, color_for_current_pixel);

			s->find_pixel_and_color(current_x, current_y,
					color_for_current_pixel, current_z, current_shading);
		}
	}
}
