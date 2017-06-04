#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>
#include <fstream>

#include "Utilities.h"
#include "TriangleOperations.h"
#include "LightingParameters.h"
#include "RenderFunctions.h"

#define MINRANGE 1
#define MAXRANGE 6
#define NUMCOLORS 8

/*double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
unsigned char RGB[NUMCOLORS][3] = {
        { 71, 71, 219 },
        { 0, 0, 91 },
        { 0, 255, 255 },
        { 0, 128, 0 },
        { 255, 255, 0 },
        { 255, 96, 0 },
        { 107, 0, 0 },
        { 224, 76, 76 } 
};*/

LightingParameters lp;

void get_color_for_vertex(double (&color)[3],float val, double *mins, double *maxs, unsigned char RGB[][3]) {
	int r;
	for (r = 0 ; r < NUMCOLORS-1 ; r++) {
		if (mins[r] <= val && val < maxs[r])
		    break;
	}

	if (r == NUMCOLORS) {
		cerr << "Could not interpolate color for " << val << endl;
		exit (EXIT_FAILURE);
	}

	double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
	color[0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
	color[1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
	color[2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
}

std::vector<Triangle> GetTriangles(const char *filename, const char *variable, const char *config) {
	int num_files, num_colors;
	double *mins;
	double *maxs;
	ifstream toRead(config);
	if (toRead.is_open()) {
		toRead >> num_files >> num_files >> num_files >> num_files >> num_files >> num_files >> num_files;
		toRead >> num_colors;
		mins = (double*)malloc((num_colors-1)*sizeof(double));
		maxs = (double*)malloc((num_colors-1)*sizeof(double));
		toRead >> mins[0];
		int i;
		for(i = 1; i < num_colors - 1; i++) {
			toRead >> mins[i];
			maxs[i-1] = mins[i];
		}
		toRead >> maxs[i-1];		
	}	
	unsigned char RGB[num_colors][3];
	int color_val;
	for(int i = 0; i < num_colors; i++) {
		for(int j = 0; j < 3; j++) {
			toRead >> color_val;
			RGB[i][j] = color_val;
		}
	}

	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName(filename);
	rdr->Update();
	if (rdr->GetOutput()->GetNumberOfCells() == 0) {
		cerr << "Unable to open file!!" << endl;
		exit (EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris     = pd->GetNumberOfCells();
	std::vector<Triangle> tris(numTris);
	vtkPoints *pts      = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
	double *color_ptr   = var->GetPointer(0);
	vtkFloatArray *n    = (vtkFloatArray *) pd->GetPointData()->GetNormals();
	float *normals      = n->GetPointer(0);
	vtkIdType npts;
	vtkIdType *ptIds;
	int idx;
        float val;
        double *pt = nullptr;
	for (idx = 0, cells->InitTraversal();
			cells->GetNextCell(npts, ptIds); idx++) {
		if (npts != 3) {
			cerr << "Non-triangles!! ???" << endl;
			exit (EXIT_FAILURE);
		}
		pt = pts->GetPoint(ptIds[0]);
		tris[idx].X[0] = pt[0];
		tris[idx].Y[0] = pt[1];
		tris[idx].Z[0] = pt[2];
		pt = pts->GetPoint(ptIds[1]);
		tris[idx].X[1] = pt[0];
		tris[idx].Y[1] = pt[1];
		tris[idx].Z[1] = pt[2];
		pt = pts->GetPoint(ptIds[2]);
		tris[idx].X[2] = pt[0];
		tris[idx].Y[2] = pt[1];
		tris[idx].Z[2] = pt[2];

		tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
		tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
		tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
		tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
		tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
		tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
		tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
		tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
		tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

		for (int j = 0; j < 3; j++) {
			val = color_ptr[ptIds[j]];
			get_color_for_vertex(tris[idx].colors[j], val, mins, maxs, RGB);
		}
	}
	return tris;
}

void transformTriangle(Triangle *t, Matrix composite, Camera camera) {
        double view_dir[3];
        double current_quadro[4];
	double transformed_vertex[4];
	for (int i = 0; i < 3; i++) {
		view_dir[0] = t->X[i] - camera.position[0];
                view_dir[1] = t->Y[i] - camera.position[1];
                view_dir[2] = t->Z[i] - camera.position[2];
		normalize_vector(view_dir);

		t->shading[i] = calculate_phong_shading(lp, view_dir, t->normals[i]);
                current_quadro[0] = t->X[i];
                current_quadro[1] = t->Y[i];
                current_quadro[2] = t->Z[i];
                current_quadro[3] = 1;
		composite.TransformPoint(current_quadro, transformed_vertex);

		if (transformed_vertex[3] != 1.) {
				transformed_vertex[0] = transformed_vertex[0] / transformed_vertex[3];
				transformed_vertex[1] = transformed_vertex[1] / transformed_vertex[3];
				transformed_vertex[2] = transformed_vertex[2] / transformed_vertex[3];
                }
		t->X[i] = transformed_vertex[0];
		t->Y[i] = transformed_vertex[1];
		t->Z[i] = transformed_vertex[2];
	}
}
