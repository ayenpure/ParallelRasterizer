#include <iostream>
#include <string.h>
#include <cmath>
#include <algorithm>
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
#include "CameraPositions.h"
#include "Utilities.h"
#include "TriangleOperations.h"
#include "MatrixOperations.h"
#include "Camera.h"
#include "Screen.h"
#include "RenderFunctions.h"

using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::pow;
using std::tan;
using std::sin;

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
	int config_id = 0,
	 		no_of_procs = 0,
			num_cameras = 0,
			read_files = 0;
	std::vector< std::vector<Triangle> > proc_triangles;
	if(argc < 4) {
		cout << "Incorrect number of arguments for execution" << endl;
		exit (EXIT_FAILURE);
	} else {
		no_of_procs = atoi(argv[1]);
		config_id = atoi(argv[2]);
		read_files = atoi(argv[3]);
		if(!read_files)
			proc_triangles = GetTrianglesForProcs(2,0,7);
	}
	double** camera_positions = get_camera_positions(config_id,&num_cameras);
	int pixels_deposited = 0;
	int pixels_deposited_per_node[no_of_procs];
	for(int file_index = 0; file_index < no_of_procs; file_index++) {
		pixels_deposited = 0;
		std::vector<Triangle> triangles;
		std::ostringstream oss;
		if(read_files) {
			oss << argv[4] << "." << file_index << ".vtk";
			triangles = GetTriangles(oss.str().c_str(), argv[4]);
			oss.str("");
			oss.clear();
		} else {
			triangles = proc_triangles[file_index];
		}
		for(int cam_index = 0; cam_index < num_cameras; cam_index++) {
			vtkImageData *image = NewImage(1000, 1000);
			unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
			int npixels = 1000 * 1000;
			for (int i = 0; i < npixels * 3; i++)
				buffer[i] = 0;
			double depth_buffer[npixels];
			for (int i = 0; i < npixels; i++)
				depth_buffer[i] = -1;
			Screen screen;
			screen.buffer = buffer;
			screen.depth_buffer = depth_buffer;
			screen.width = 1000;
			screen.height = 1000;
			double camera_position[3], focus_point[3];
			get_camera_position_and_focus(config_id, cam_index, num_cameras, camera_positions,
				camera_position, focus_point);
			Camera camera = GetCamera(camera_position, focus_point);
			Matrix camera_transform = camera.CameraTransform();
			Matrix view_transform = camera.ViewTransform();
			Matrix device_transform = camera.DeviceTransform(screen);
			Matrix composite = get_total_transform_matrix(camera_transform,
					view_transform, device_transform);
			for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
				Triangle t = triangles[vecIndex];
				transformTriangle(&t, composite, camera);
				if (t.is_flat_bottom_triangle()) {
					pixels_deposited += scan_line(&t, &screen);
				} else {
					Triangle t1, t2;
					t.split_triangle(&t1, &t2);
					pixels_deposited += scan_line(&t1, &screen);
					pixels_deposited += scan_line(&t2, &screen);
				}
			}
			oss << "camera_" << file_index << "_" << cam_index;
			WriteImage(image, oss.str().c_str());
			oss.str("");
			oss.clear();
			free(buffer);
		}
		pixels_deposited_per_node[file_index] = pixels_deposited / num_cameras;
	}
	free(camera_positions);
	cout << "\n\n";
	for(int i = 0; i < no_of_procs; i++) {
		cout << pixels_deposited_per_node[i] << endl;
	}
}
