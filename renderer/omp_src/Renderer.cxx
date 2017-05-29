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
#include <omp.h>

#define HEIGHT 1000
#define WIDTH 1000

using std::cerr;
using std::endl;

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
	if(argc < 2) {
		cout << "Usage : render <mesh file> <variable/field>" << endl;
		exit (EXIT_FAILURE);
	} else {
		triangles = GetTriangles(argv[1], argv[2]);
	}

	int no_of_triangles = triangles.size();

	vtkImageData *image = NewImage(HEIGHT, WIDTH);
	unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
	int npixels = WIDTH * HEIGHT;

        #pragma omp parallel for num_threads(4)
	for (int i = 0; i < npixels * 3; i++)
		buffer[i] = 0;

	double *depth_buffer = (double*)malloc(npixels*sizeof(double));

        #pragma omp parallel for num_threads(4)
	for (int i = 0; i < npixels; i++)
		depth_buffer[i] = -1;

	Screen screen;
	screen.buffer = buffer;
	screen.depth_buffer = depth_buffer;
	screen.width = WIDTH;
	screen.height = HEIGHT;
	double camera_position[] = {0,40,40};
	double focus_point[] = {0,0,0};
		
	Camera camera = GetCamera(camera_position, focus_point);
		
	Matrix camera_transform = camera.CameraTransform();
	Matrix view_transform = camera.ViewTransform();
	Matrix device_transform = camera.DeviceTransform(screen);
	Matrix composite = get_total_transform_matrix(camera_transform,
			view_transform, device_transform);
	
        #pragma omp parallel for num_threads(3)
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

	std::ostringstream oss;
	oss << "outputimage";
	WriteImage(image, oss.str().c_str());
	oss.clear();
}