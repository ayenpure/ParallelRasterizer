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
NewRGBImage(int width, int height) {
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	return img;
}

vtkImageData *
NewRGBAImage(int width, int height) {
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
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
	int height = 1024, width = 1024;
	int config_id = 0,
	 		no_of_procs = 0,
			num_cameras = 0,
			read_files = 0;
	std::vector<Triangle> triangles;
	if(argc < 4) {
		cout << "Incorrect number of arguments for execution\n"
		<< "Usage : <program> <no of pros> <camera config> <read flag>" << endl;
		exit (EXIT_FAILURE);
	} else {
		no_of_procs = atoi(argv[1]);
		config_id = atoi(argv[2]);
		read_files = atoi(argv[3]);
		if(!read_files)
			triangles = GetTriangles(2,0,7);
		else {
			triangles = GetTrianglesFromFiles(no_of_procs, argv[4]);
		}
	}
	//exit (EXIT_FAILURE);
	double** camera_positions = get_camera_positions(config_id,&num_cameras);
	int active_pixels[num_cameras];
	int no_of_triangles = triangles.size();
	for (int i = 0; i < num_cameras; i++)
		active_pixels[i] = 0;
	for (int cam_index = 0; cam_index < num_cameras; cam_index++) {
		vtkImageData *image;
		unsigned char *buffer;
		image = NewRGBImage(height, width);
		int npixels = height * width;
		buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
		double *depth_buffer = (double*)malloc(npixels*sizeof(double));
		for (int i = 0; i < npixels * 3; i++)
			buffer[i] = 0;
		for (int i = 0; i < npixels; i++)
			depth_buffer[i] = -2;
		Screen screen;
		screen.buffer = buffer;
		screen.depth_buffer = depth_buffer;
		screen.width = width;
		screen.height = height;
		double camera_position[3], focus_point[3];
		get_camera_position_and_focus(config_id, cam_index, num_cameras,
			 camera_positions, camera_position, focus_point);
		Camera camera = GetCamera(camera_position, focus_point);
		Matrix camera_transform = camera.CameraTransform();
		//cout<<"\nCamera Transform Matrix :\n";camera_transform.Print(std::cout);
		Matrix view_transform = camera.ViewTransform();
		//cout<<"\nView Transform Matrix :\n";view_transform.Print(std::cout);
		Matrix device_transform = camera.DeviceTransform(screen);
		//cout<<"\nDevice Transform Matrix :\n";device_transform.Print(std::cout);
		Matrix composite = get_total_transform_matrix(camera_transform,
				view_transform, device_transform);
		//cout<<"\nComposite Matrix :\n";composite.Print(std::cout);
		for (int vecIndex = 0; vecIndex < triangles.size(); vecIndex++) {
			Triangle t = triangles[vecIndex];
			transformTriangle(&t, composite, camera);
			if (t.is_flat_bottom_triangle()) {
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t, &screen);
			} else {
				Triangle t1, t2;
				t.split_triangle(&t1, &t2);
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t1, &screen);
				/*triangle_pixels[vecIndex] = triangle_pixels[vecIndex] + */scan_line(&t2, &screen);
			}
		}
		std::ostringstream oss;
		oss << "camera_" << cam_index;
		active_pixels[cam_index] = 0;
		for(int pxl =0; pxl < 3*npixels; pxl+=3) {
			if(!(buffer[pxl] == 0 &&
				buffer[pxl+1] == 0 &&
				buffer[pxl+2] == 0 ))
				active_pixels[cam_index]++;
		}
		active_pixels[cam_index] = active_pixels[cam_index] / no_of_procs;
		WriteImage(image, oss.str().c_str());
		oss.str("");
		oss.clear();
		free(buffer);
		image = NewRGBAImage(height,width);
		buffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
		/*For Shawn*/
		for (int i = 0; i < npixels * 4; i++)
			buffer[i] = 0;
		for(int i = 0; i < npixels; i ++) {
			if(!(depth_buffer[i] >= -1 && depth_buffer[i] <= 1))
				continue;
			double value = depth_buffer[i];
			double newValue = range_transform(0.,1.,1.,-1., value);
			double rgba[] = {0,0,0,0};
			getRGBAforDepth(newValue, rgba);
			buffer[i*4] = rgba[0]*255;
			buffer[i*4 + 1] = rgba[1]*255;
			buffer[i*4 + 2] = rgba[2]*255;
			buffer[i*4 + 3] = rgba[3]*255;
		}
		oss << "depth_" << cam_index;
		WriteImage(image, oss.str().c_str());
		oss.str("");
		oss.clear();
		free(buffer);
		free(depth_buffer);
	}
	for (int i = 0; i < num_cameras; i++) {
		cout << active_pixels[i] << endl;
	}
}
