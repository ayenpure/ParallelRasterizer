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

#define MINRANGE 0.0653
#define MAXRANGE 98.52
#define NUMCOLORS 5

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

struct vert {
  double x,y,z;
  bool operator==(const vert &other) const {
    return ( x==other.x &&
             y==other.y &&
             z==other.z );
  }
};

namespace std {
  template<> struct hash<vert>
  {
    /*typedef vertex argument_type;
    typedef std::size_t result_type;*/
    size_t operator()(vert const& v) const {
      size_t const h1 ( std::hash<double>{}(v.x) );
      size_t const h2 ( std::hash<double>{}(v.y) );
      size_t const h3 ( std::hash<double>{}(v.z) );
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };
}

/*void process_for_vertex_normals(std::vector<Triangle> tris) {
	unordered_map<vert, pair<double*, int>> vertices;
	int no_of_triangles = tris.size();
  std::cout << "number of triangels : " << no_of_triangles << endl;
	for(int i = 0; i <= no_of_triangles; i++) {
		double normal[3];
		bool normal_calculated = tris[i].calculate_normal(normal);
		if(!normal_calculated) {
			std::cout << "Normal not calculated properly at index " << i << endl;
			continue;
		}
		for(int j = 0; j < 3; j++) {
			//Create entry in vertices for each vertex
			vert v = {tris[i].X[j],tris[i].Y[j],tris[i].Z[j]};

			int count = vertices.count(v);
			if(count == 0) {
        // calculate triangle normal
        // add vertex to the map as key and the normal as value
        // count in the pair is 1
				double* temp = (double*)malloc(3*sizeof(double));
				vector_copy(normal, temp);
        pair<double*, int> value = std::make_pair(temp, 1);
        pair<vert, pair<double*, int>> toInsert = std::make_pair(v, value);
        vertices.insert(toInsert);
			} else {
				// fetch the existing data
        // calculate triangle normal
				// add the normal to existing normal
        // increase count in the pair
        // push new value (as it's a vector, it would change in actual values)
        // doubtful about the count though
        auto toProcess = vertices.find(v);
        if(toProcess != vertices.end()) {
					double* temp = (double*)malloc(3*sizeof(double));
          auto value = toProcess->second;
          temp[0] = value.first[0] + normal[0];
          temp[1] = value.first[1] + normal[1];
          temp[2] = value.first[2] + normal[2];
          int count = ++value.second;
					free(value.first);
					pair<double*, int> newValue = std::make_pair(temp, count);
					toProcess->second = newValue;
        }
			}
		}
	}
  std::cout << "number of vertices : " << vertices.size() << endl;
	std::cout << "mymap contains : " << endl;
	for ( auto toPrint = vertices.begin(); toPrint != vertices.end(); ++toPrint ) {
		 auto addedNormals = toPrint->second;
		 double* normals = addedNormals.first;
		 int count = addedNormals.second;
		 normals[0] /= count;
		 normals[1] /= count;
		 normals[2] /= count;
		 //normalize_vector(normals);
	}
	for(int i = 0 ;i < no_of_triangles; i++) {
		for(int j = 0; j < 3; j++) {
			vert v = {tris[i].X[j],tris[i].Y[j],tris[i].Z[j]};
			auto toProcess = vertices.find(v);
			if(toProcess != vertices.end()) {
				tris[i].normals[j][0] = toProcess->second.first[0];
				tris[i].normals[j][1] = toProcess->second.first[1];
				tris[i].normals[j][2] = toProcess->second.first[2];
			}
		}
	}
}*/

std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
{
		std::vector<Triangle> output(4*list.size());
		for (unsigned int i = 0 ; i < list.size() ; i++)
		{
				double mid1[3], mid2[3], mid3[3];
				mid1[0] = (list[i].X[0]+list[i].X[1])/2;
				mid1[1] = (list[i].Y[0]+list[i].Y[1])/2;
				mid1[2] = (list[i].Z[0]+list[i].Z[1])/2;
				mid2[0] = (list[i].X[1]+list[i].X[2])/2;
				mid2[1] = (list[i].Y[1]+list[i].Y[2])/2;
				mid2[2] = (list[i].Z[1]+list[i].Z[2])/2;
				mid3[0] = (list[i].X[0]+list[i].X[2])/2;
				mid3[1] = (list[i].Y[0]+list[i].Y[2])/2;
				mid3[2] = (list[i].Z[0]+list[i].Z[2])/2;

				output[4*i+0].X[0] = list[i].X[0];
				output[4*i+0].Y[0] = list[i].Y[0];
				output[4*i+0].Z[0] = list[i].Z[0];
				output[4*i+0].X[1] = mid1[0];
				output[4*i+0].Y[1] = mid1[1];
				output[4*i+0].Z[1] = mid1[2];
				output[4*i+0].X[2] = mid3[0];
				output[4*i+0].Y[2] = mid3[1];
				output[4*i+0].Z[2] = mid3[2];

				for (int j = 0; j < 3; j++) {
					output[4*i+0].colors[j][0] = 0./255.0;
					output[4*i+0].colors[j][1] = 96./255.0;
					output[4*i+0].colors[j][2] = 69./255.0;
				}


				output[4*i+1].X[0] = list[i].X[1];
				output[4*i+1].Y[0] = list[i].Y[1];
				output[4*i+1].Z[0] = list[i].Z[1];
				output[4*i+1].X[1] = mid2[0];
				output[4*i+1].Y[1] = mid2[1];
				output[4*i+1].Z[1] = mid2[2];
				output[4*i+1].X[2] = mid1[0];
				output[4*i+1].Y[2] = mid1[1];
				output[4*i+1].Z[2] = mid1[2];

				for (int j = 0; j < 3; j++) {
					output[4*i+1].colors[j][0] = 0./255.0;
					output[4*i+1].colors[j][1] = 96./255.0;
					output[4*i+1].colors[j][2] = 69./255.0;
				}

				output[4*i+2].X[0] = list[i].X[2];
				output[4*i+2].Y[0] = list[i].Y[2];
				output[4*i+2].Z[0] = list[i].Z[2];
				output[4*i+2].X[1] = mid3[0];
				output[4*i+2].Y[1] = mid3[1];
				output[4*i+2].Z[1] = mid3[2];
				output[4*i+2].X[2] = mid2[0];
				output[4*i+2].Y[2] = mid2[1];
				output[4*i+2].Z[2] = mid2[2];

				for (int j = 0; j < 3; j++) {
					output[4*i+2].colors[j][0] = 0./255.0;
					output[4*i+2].colors[j][1] = 96./255.0;
					output[4*i+2].colors[j][2] = 69./255.0;
				}

				output[4*i+3].X[0] = mid1[0];
				output[4*i+3].Y[0] = mid1[1];
				output[4*i+3].Z[0] = mid1[2];
				output[4*i+3].X[1] = mid2[0];
				output[4*i+3].Y[1] = mid2[1];
				output[4*i+3].Z[1] = mid2[2];
				output[4*i+3].X[2] = mid3[0];
				output[4*i+3].Y[2] = mid3[1];
				output[4*i+3].Z[2] = mid3[2];

				for (int j = 0; j < 3; j++) {
					output[4*i+3].colors[j][0] = 0./255.0;
					output[4*i+3].colors[j][1] = 96./255.0;
					output[4*i+3].colors[j][2] = 69./255.0;
				}
		}
		return output;
}

std::vector<Triangle> SplitTriangle(std::vector<Triangle> list, int parts) {
	if(parts == 0)
		return list;
	std::vector<Triangle> octant_triangles(list.size()*parts);
	int count = 0;
	for(int i = 0; i < list.size(); i++) {
		Triangle t = list[i];
		int offset_index = (t.Y[0] > t.Y[1] && t.Y[0] > t.Y[2]) ? 0 :
					(t.Y[1] > t.Y[0] && t.Y[1] > t.Y[2]) ? 1 : 2;
		int adj_1, adj_2;
		adj_1 = (offset_index + 1) % 3;
		adj_2 = (offset_index + 2) % 3;
		for(int j = 0; j < parts; j++) {
			Triangle tpart;
			double point[3];
			point[0] = interpolate(0, parts, t.X[adj_1], t.X[adj_2], j+1);
			point[1] = interpolate(0, parts, t.Y[adj_1], t.Y[adj_2], j+1);
			point[2] = interpolate(0, parts, t.Z[adj_1], t.Z[adj_2], j+1);
			if(j == 0) {
				tpart.X[0] = t.X[offset_index];
				tpart.X[1] = t.X[adj_1];
				tpart.X[2] = point[0];
				tpart.Y[0] = t.Y[offset_index];
				tpart.Y[1] = t.Y[adj_1];
				tpart.Y[2] = point[1];
				tpart.Z[0] = t.Z[offset_index];
				tpart.Z[1] = t.Z[adj_1];
				tpart.Z[2] = point[2];
			} else {
				tpart.X[0] = t.X[offset_index];
				tpart.X[1] = octant_triangles[count-1].X[2];
				tpart.X[2] = point[0];
				tpart.Y[0] = t.Y[offset_index];
				tpart.Y[1] = octant_triangles[count-1].Y[2];
				tpart.Y[2] = point[1];
				tpart.Z[0] = t.Z[offset_index];
				tpart.Z[1] = octant_triangles[count-1].Z[2];
				tpart.Z[2] = point[2];
			}
			octant_triangles[count++] = tpart;
		}
	}
	cout << "number of octant triangles : " << count << endl;
	return octant_triangles;
}

Triangle rotate_triangle(Triangle t, double angle, char axis) {
	Triangle temp;
	for (int i = 0; i < 3; i++) {
		double current_vertex[3] = { t.X[i], t.Y[i], t.Z[i] };
		double transformed_vertex[3] = {0, 0, 0};
		rotate(current_vertex, angle, axis,transformed_vertex);
		temp.X[i] = transformed_vertex[0];
		temp.Y[i] = transformed_vertex[1];
		temp.Z[i] = transformed_vertex[2];
	}
	return temp;
}

std::vector<Triangle> get_all_octants(Triangle t) {
	int num_octants = 8;
	std::vector<Triangle> octants(num_octants);
	octants[0] = t;
	for(int i = 1; i < num_octants; i++) {
		if(i == 4 )
			octants[i] = rotate_triangle(octants[0],M_PI / 2, 'z');
		else
			octants[i] = rotate_triangle(octants[i-1],M_PI / 2, 'y');
	}
	return octants;
}

std::vector<Triangle> GetTriangles(int split_parts, int split_rec, int tri_grain) {
	std::vector<Triangle> list;
	Triangle t;
	t.X[0] = 1;
	t.Y[0] = 0;
	t.Z[0] = 0;
	t.X[1] = 0;
	t.Y[1] = 1;
	t.Z[1] = 0;
	t.X[2] = 0;
	t.Y[2] = 0;
	t.Z[2] = 1;
	for (int j = 0; j < 3; j++) {
		t.colors[j][0] = 0./255.0;
		t.colors[j][1] = 96./255.0;
		t.colors[j][2] = 69./255.0;
	}
	list = get_all_octants(t);
	list = SplitTriangle(list, split_parts);
	for(int i = 0; i < split_rec; i++)
		list = SplitTriangle(list);
	cout << "To split for processor work : " << list.size() << endl;
	std::vector< std::vector<Triangle> > proc_parted_triangles(list.size());
	for(int i = 0; i < list.size(); i++) {
		std::vector<Triangle> toExpand;
		toExpand.push_back(list[i]);
		for(int j = 0; j < tri_grain; j++) {
			toExpand = SplitTriangle(toExpand);
		}
		for(int j = 0; j < toExpand.size(); j++) {
			Triangle t = toExpand[j];
			for(int k = 0; k < 3; k++) {
				double ptMag = sqrt(t.X[k]*t.X[k]+
														t.Y[k]*t.Y[k]+
														t.Z[k]*t.Z[k]);
				t.X[k] = (t.X[k] / ptMag)*10;
				t.Y[k] = (t.Y[k] / ptMag)*10;
				t.Z[k] = (t.Z[k] / ptMag)*10;
			}
			for (int k = 0; k < 3; k++) {
				t.colors[k][0] = 0./255.0;
				t.colors[k][1] = 96./255.0;
				t.colors[k][2] = 69./255.0;
			}
			t.calculate_normals();
			toExpand[j] = t;
		}
		proc_parted_triangles[i] = toExpand;
	}
	int list_out = (split_rec == 0) ? 8 * split_parts
																			: 8 * split_parts * pow(4, split_rec);
	int total_triangles = (tri_grain == 0) ? list_out : list_out * pow(4, tri_grain);
	std::vector<Triangle> newlist(total_triangles);
	for(int j = 0; j < list_out; j++) {
		newlist.insert(newlist.end(), proc_parted_triangles[j].begin(), proc_parted_triangles[j].end());
	}
	cout << "Triangles to process : " << newlist.size() << endl;
	return newlist;
}

std::vector< std::vector<Triangle> > GetTrianglesForProcs(int split_parts, int split_rec, int tri_grain) {
	std::vector<Triangle> list;
	Triangle t;
	t.X[0] = 1;
	t.Y[0] = 0;
	t.Z[0] = 0;
	t.X[1] = 0;
	t.Y[1] = 1;
	t.Z[1] = 0;
	t.X[2] = 0;
	t.Y[2] = 0;
	t.Z[2] = 1;
	for (int j = 0; j < 3; j++) {
		t.colors[j][0] = 0./255.0;
		t.colors[j][1] = 96./255.0;
		t.colors[j][2] = 69./255.0;
	}
	list = get_all_octants(t);
	list = SplitTriangle(list, split_parts);
	for(int i = 0; i < split_rec; i++)
		list = SplitTriangle(list);
	cout << "To split for processor work : " << list.size() << endl;
	std::vector< std::vector<Triangle> > proc_parted_triangles(list.size());
	for(int i = 0; i < list.size(); i++) {
		std::vector<Triangle> toExpand;
		toExpand.push_back(list[i]);
		for(int j = 0; j < tri_grain; j++) {
			toExpand = SplitTriangle(toExpand);
		}
		for(int j = 0; j < toExpand.size(); j++) {
			Triangle t = toExpand[j];
			for(int k = 0; k < 3; k++) {
				double ptMag = sqrt(t.X[k]*t.X[k]+
														t.Y[k]*t.Y[k]+
														t.Z[k]*t.Z[k]);
				t.X[k] = (t.X[k] / ptMag)*10;
				t.Y[k] = (t.Y[k] / ptMag)*10;
				t.Z[k] = (t.Z[k] / ptMag)*10;
			}
			for (int k = 0; k < 3; k++) {
				t.colors[k][0] = 0./255.0;
				t.colors[k][1] = 96./255.0;
				t.colors[k][2] = 69./255.0;
			}
			t.calculate_normals();
			toExpand[j] = t;
		}
		proc_parted_triangles[i] = toExpand;
		//cout << "Proc " << i << " has " << toExpand.size() << " triangles" << endl;
	}
	/*std::vector<Triangle> newlist(64*1024);
	for(int j = 0; j < 64; j++) {
		newlist.insert(newlist.end(), proc_parted_triangles[j].begin(), proc_parted_triangles[j].end());
	}
	cout << "Triangles to process : " << newlist.size() << endl;
	return newlist;*/
	return proc_parted_triangles;
}

void get_color_for_vertex(double* color,float val) {
	int num_colors = NUMCOLORS;
  double range_inc = (MAXRANGE - MINRANGE) / num_colors;
	/*double mins[num_colors-1] = { 1, 2.25, 3.5, 4.75};
	double maxs[num_colors-1] = { 2.25, 3.5, 4.75, 6};*/
	double mins[num_colors-1] = { MINRANGE,MINRANGE+range_inc,MINRANGE+2*range_inc,MINRANGE+3*range_inc};
	double maxs[num_colors-1] = { MINRANGE+range_inc,MINRANGE+2*range_inc,MINRANGE+3*range_inc, MAXRANGE};
	unsigned char RGB[num_colors][3] = {
		{0, 0, 255},
		{0, 204, 255},
		{0, 153, 0},
		{255, 204, 0},
		{255, 0, 0},
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

std::vector<Triangle> GetTriangles(const char *filename, char *variable) {
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
	//vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
	//double *color_ptr = var->GetPointer(0);
	vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray(variable);
	float *color_ptr = var->GetPointer(0);
	vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
	//float *normals = n->GetPointer(0);
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

		for (int j = 0; j < 3; j++) {
			float val = color_ptr[ptIds[j]];
			get_color_for_vertex(tris[idx].colors[j], val);
		}
		tris[idx].calculate_normals();
	}
	return tris;
}

std::vector<Triangle> GetTrianglesFromFiles(int no_of_procs, char *variable) {
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
			/*tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
			 tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
			 tris[idx].normals[0][2] = normals[3*ptIds[0]+2];*/
			pt = pts->GetPoint(ptIds[1]);
			tris[idx].X[1] = pt[0];
			tris[idx].Y[1] = pt[1];
			tris[idx].Z[1] = pt[2];
			/*tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
			 tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
			 tris[idx].normals[1][2] = normals[3*ptIds[1]+2];*/
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
}

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
