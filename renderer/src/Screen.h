#ifndef SCREEN_H_
#define SCREEN_H_

#include "Utilities.h"

using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

extern int numFrags;
class Screen {
public:
	unsigned char *buffer;
	double *depth_buffer;
	int width, height;

	bool find_pixel_and_color(int x, int y, double *color, double current_depth,
			double shading_amount) {
		/*
                numFrags++;
                if (numFrags % 10000 == 0)
                      cerr << "Num frags = " << numFrags << endl;
		*/
		/*
		 * Ensure the pixels to be painted are in the frame.
		 */
		if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
			int buffer_index = (y * 3 * width) + (x * 3);
			int depth_buffer_index = y * width + x;
			if (buffer_index < width * height * 3
					&& (current_depth >= -1 && current_depth <= 1)
					&& current_depth >= depth_buffer[depth_buffer_index]) {
				/*double shading_amount = calculate_phong_shading(lp,
				 view_direction, current_normal);*/
				buffer[buffer_index++] = min(
						ceil441((shading_amount * color[0]) * 255),
						(double) 255);
				buffer[buffer_index++] = min(
						ceil441((shading_amount * color[1]) * 255),
						(double) 255);
				buffer[buffer_index] = min(
						ceil441((shading_amount * color[2]) * 255),
						(double) 255);
				depth_buffer[depth_buffer_index] = current_depth;
				return true;
			}
		}
		return false;
	}

	/*
	 * This method calculates the color for a certain pixel in a trangle given
	 * its intercepts for a scanline and the colors at the intercepts.
	 */
	void calculate_color_for_pixel(double left_intercept,
			double right_intercept, double current_x_pos,
			double *color_at_left_intercept, double *color_at_right_intercept,
			double *current_color) {

		current_color[0] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[0], color_at_right_intercept[0],
				current_x_pos);
		current_color[1] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[1], color_at_right_intercept[1],
				current_x_pos);
		current_color[2] = interpolate(left_intercept, right_intercept,
				color_at_left_intercept[2], color_at_right_intercept[2],
				current_x_pos);
	}
};

#endif
