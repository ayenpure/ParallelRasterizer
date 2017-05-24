#ifndef RENDERFUNCTIONS_H_
#define RENDERFUNCTIONS_H_

#include "LightingParameters.h"
#include "TriangleOperations.h"

double calculate_for_specular_lighting(LightingParameters lp, double *V,
		double* N);

double calculate_phong_shading(LightingParameters lp, double *view_direction,
		double *normal);

double scan_line(Triangle *t, Screen *s);

#endif
