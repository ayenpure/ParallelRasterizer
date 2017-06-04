#include <iostream>
#include <string.h>
#include <vector>
#include <sstream>
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
#include "MatrixOperations.h"
#include "Camera.h"
#include "Screen.h"
#include "RenderFunctions.h"
#include <fstream>
#include <ctime>
#include <ratio>
#include <chrono>

#define HEIGHT 1000
#define WIDTH 1000

using namespace std;
using namespace std::chrono;

vtkImageData *
NewImage(int width, int height) {
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	return img;
}

void WriteImage(vtkImageData *img, const char *filename) {
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

int main(int argc, char *argv[]) {
	std::vector<Triangle> triangles;
	if(argc < 3) {
		cout << "Usage : render <mesh file> <variable/field> <config>" << endl;
		exit (EXIT_FAILURE);
	} else {
		triangles = GetTriangles(argv[1], argv[2], argv[3]);
	}

	double camera_position[3];
	double focus_point[3];
	ifstream toRead(argv[3]);
	if (toRead.is_open()) {
		toRead >> camera_position[0] >> camera_position[1] >> camera_position[2];
		toRead >> focus_point[0] >> focus_point[1] >> focus_point[3];
		toRead.close();
	}

	int no_of_triangles = triangles.size();

	vtkImageData *image = NewImage(HEIGHT, WIDTH);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = WIDTH * HEIGHT;
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;
	double *depth_buffer = (double*)malloc(npixels*sizeof(double));
	for (int i = 0; i < npixels; i++)
		depth_buffer[i] = -1;
	Screen screen;
	screen.buffer = buffer;
	screen.depth_buffer = depth_buffer;
	screen.width = WIDTH;
	screen.height = HEIGHT;
		
	Camera camera = GetCamera(camera_position, focus_point);
		
	Matrix camera_transform = camera.CameraTransform();
	Matrix view_transform = camera.ViewTransform();
	Matrix device_transform = camera.DeviceTransform(screen);
	Matrix composite = get_total_transform_matrix(camera_transform,
			view_transform, device_transform);

	high_resolution_clock::time_point r_start = high_resolution_clock::now();

	for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
				
		Triangle t = triangles[vecIndex];

		transformTriangle(&t, composite, camera);
		
		if (t.is_flat_bottom_triangle()) {
			scan_line(&t, &screen);
		} else {
			Triangle t1, t2;
			t.split_triangle(&t1, &t2);
			scan_line(&t1, &screen);
			scan_line(&t2, &screen);
		}
	}
	
	high_resolution_clock::time_point r_end = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(r_end - r_start);
	std::cout << "Time to render image : " << time_span.count() << " seconds" << endl;

	std::ostringstream oss;
	oss << "outputimage";
	WriteImage(image, oss.str().c_str());
	oss.clear();
}
