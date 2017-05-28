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

#include "Utilities.h"
#include "TriangleOperations.h"
#include "LightingParameters.h"
#include "RenderFunctions.h"

#define MINRANGE 1
#define MAXRANGE 6
#define NUMCOLORS 8

using std::cout;
using std::unordered_map;
using std::pair;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

LightingParameters lp;

void get_color_for_vertex(double* color,float val) {
	int num_colors = NUMCOLORS;
        //double range_inc = (MAXRANGE - MINRANGE) / num_colors;
	/*double mins[num_colors-1] = { 1, 2.25, 3.5, 4.75};
	double maxs[num_colors-1] = { 2.25, 3.5, 4.75, 6};*/
	//double mins[num_colors-1] = { MINRANGE,MINRANGE+range_inc,MINRANGE+2*range_inc,MINRANGE+3*range_inc};
	//double maxs[num_colors-1] = { MINRANGE+range_inc,MINRANGE+2*range_inc,MINRANGE+3*range_inc, MAXRANGE};
	double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
	double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
	unsigned char RGB[num_colors][3] = {
		{ 71, 71, 219 },
		{ 0, 0, 91 },
		{ 0, 255, 255 },
		{ 0, 128, 0 },
		{ 255, 255, 0 },
		{ 255, 96, 0 },
		{ 107, 0, 0 },
		{ 224, 76, 76 } 
	};
	int r;
	for (r = 0 ; r < num_colors-1 ; r++) {
		if (mins[r] <= val && val < maxs[r])
		break;
	}
	if (r == num_colors) {
		cerr << "Could not interpolate color for " << val << endl;
		exit (EXIT_FAILURE);
	}
	double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
	color[0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
	color[1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
	color[2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
	/*color[0] = 1;
	color[1] = 0;
	color[2] = 0;*/
}

std::vector<Triangle> GetTriangles(const char *filename,const char *variable) {
	vtkPolyDataReader *rdr = vtkPolyDataReader::New();
	rdr->SetFileName(filename);
	//cerr << "Reading" << endl;
	rdr->Update();
	//cerr << "Done reading" << endl;
	if (rdr->GetOutput()->GetNumberOfCells() == 0) {
		cerr << "Unable to open file!!" << endl;
		exit (EXIT_FAILURE);
	}
	vtkPolyData *pd = rdr->GetOutput();
	int numTris = pd->GetNumberOfCells();
	std::vector<Triangle> tris(numTris);
	vtkPoints *pts = pd->GetPoints();
	vtkCellArray *cells = pd->GetPolys();
	vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
	double *color_ptr = var->GetPointer(0);
	//vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray(variable);
	//float *color_ptr = var->GetPointer(0);
	vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
	float *normals = n->GetPointer(0);
	vtkIdType npts;
	vtkIdType *ptIds;
	int idx;
	for (idx = 0, cells->InitTraversal();
			cells->GetNextCell(npts, ptIds); idx++) {
		if (npts != 3) {
			cerr << "Non-triangles!! ???" << endl;
			exit (EXIT_FAILURE);
		}
		double *pt = NULL;
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
			float val = color_ptr[ptIds[j]];
			get_color_for_vertex(tris[idx].colors[j], val);
		}
		//tris[idx].calculate_normals();
	}
	return tris;
}

/*std::vector<Triangle> GetTrianglesFromFiles(int no_of_procs, char *variable) {
	int index = 0;
	std::vector<Triangle> tris(0);
	for (int file_index = 0; file_index < no_of_procs; file_index++) {
		vtkPolyDataReader *rdr = vtkPolyDataReader::New();
		std::ostringstream oss;
		oss << variable << "." << file_index << ".vtk";
		rdr->SetFileName(oss.str().c_str());
		oss.str("");
		oss.clear();
		//cerr << "Reading :" << file_index << endl;
		rdr->Update();
		//cerr << "Done reading" << endl;
		if (rdr->GetOutput()->GetNumberOfCells() == 0) {
			cerr << "Unable to open file!!" << endl;
			exit (EXIT_FAILURE);
		}
		vtkPolyData *pd = rdr->GetOutput();
		int numTris = pd->GetNumberOfCells();
		tris.resize(tris.size() + numTris);
		vtkPoints *pts = pd->GetPoints();
		vtkCellArray *cells = pd->GetPolys();
		//vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
		//double *color_ptr = var->GetPointer(0);
		vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray(variable);
		float *color_ptr = var->GetPointer(0);
		vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
		//float *normals = n->GetPointer(0);
		vtkIdType npts;
		vtkIdType *ptIds;
		int idx;
		for (idx = index, cells->InitTraversal();
				cells->GetNextCell(npts, ptIds); idx++, index++) {
			if (npts != 3) {
				cerr << "Non-triangles!! ???" << endl;
				exit (EXIT_FAILURE);
			}
			double *pt = NULL;
			pt = pts->GetPoint(ptIds[0]);
			tris[idx].X[0] = pt[0];
			tris[idx].Y[0] = pt[1];
			tris[idx].Z[0] = pt[2];
			tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
			 tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
			 tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
			pt = pts->GetPoint(ptIds[1]);
			tris[idx].X[1] = pt[0];
			tris[idx].Y[1] = pt[1];
			tris[idx].Z[1] = pt[2];
			tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
			 tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
			 tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
			pt = pts->GetPoint(ptIds[2]);
			tris[idx].X[2] = pt[0];
			tris[idx].Y[2] = pt[1];
			tris[idx].Z[2] = pt[2];
			for (int j = 0; j < 3; j++) {
				float val = color_ptr[ptIds[j]];
				get_color_for_vertex(tris[idx].colors[j], val);
			}
			tris[idx].calculate_normals();
		}
	}
	//process_for_vertex_normals(tris);
	return tris;
}*/

void transformTriangle(Triangle *t, Matrix composite, Camera camera) {
	for (int i = 0; i < 3; i++) {
		double view_dir[3] = { t->X[i] - camera.position[0], t->Y[i]
				- camera.position[1], t->Z[i] - camera.position[2] };
		normalize_vector(view_dir);
		t->shading[i] = calculate_phong_shading(lp, view_dir, t->normals[i]);
		double current_quadro[4] = { t->X[i], t->Y[i], t->Z[i], 1. };
		double transformed_vertex[4];
		composite.TransformPoint(current_quadro, transformed_vertex);
		if (transformed_vertex[3] != 1.) {
			for (int j = 0; j < 3; j++)
				transformed_vertex[j] = transformed_vertex[j] / transformed_vertex[3];
		}
		t->X[i] = transformed_vertex[0];
		t->Y[i] = transformed_vertex[1];
		t->Z[i] = transformed_vertex[2];
	}
}
