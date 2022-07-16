// Copyright (c) 2022 ETH Zurich, Christian R. Steger
// MIT License

#include "shadow_comp.h"
#include <cstdio>
#include <embree3/rtcore.h>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <stdio.h>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <string.h>
#include <tbb/parallel_for.h>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace shapes;

//#############################################################################
// Auxiliary functions
//#############################################################################

// Compute linear index from subscripts (2D-array)
inline size_t lin_ind_2d(size_t dim_1, size_t ind_0, size_t ind_1) {
	return (ind_0 * dim_1 + ind_1);
}

// Compute linear index from subscripts (3D-array)
inline size_t lin_ind_3d(size_t dim_1, size_t dim_2, size_t ind_0, size_t ind_1,
	size_t ind_2) {
	return (ind_0 * (dim_1 * dim_2) + ind_1 * dim_2 + ind_2);
}

// Convert degree to radian
inline float deg2rad(float ang) {
	return ((ang / 180.0) * M_PI);
}

// Convert radian to degree
inline float rad2deg(float ang) {
	return ((ang / M_PI) * 180.0);
}

// Cross product
inline void cross_prod(float a_x, float a_y, float a_z, float b_x, float b_y,
	float b_z, float &c_x, float &c_y, float &c_z) {
	c_x = a_y * b_z - a_z * b_y;
    c_y = a_z * b_x - a_x * b_z;
    c_z = a_x * b_y - a_y * b_x;
}

// Matrix vector multiplication
inline void mat_vec_mult(float (&mat)[3][3], float (&vec)[3],
	float (&vec_res)[3]) {

	vec_res[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2];
    vec_res[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2];
    vec_res[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2];

}

//#############################################################################
// Miscellaneous
//#############################################################################

// Namespace
#if defined(RTC_NAMESPACE_USE)
	RTC_NAMESPACE_USE
#endif

// Error function
void errorFunction(void* userPtr, enum RTCError error, const char* str) {
	printf("error %d: %s\n", error, str);
}

// Initialisation of device and registration of error handler
RTCDevice initializeDevice() {
	RTCDevice device = rtcNewDevice(NULL);
  	if (!device) {
    	printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));
    }
  	rtcSetDeviceErrorFunction(device, errorFunction, NULL);
  	return device;
}

//#############################################################################
// Create scene from geometries
//#############################################################################

// Structures for triangle and quad
struct Triangle { int v0, v1, v2; };
struct Quad { int v0, v1, v2, v3; };
// -> above structures must contain 32-bit integers (-> Embree documentation).
//    Theoretically, these integers should be unsigned but the binary
//    representation until 2'147'483'647 is identical between signed/unsigned
//    integer.

// Initialise scene
RTCScene initializeScene(RTCDevice device, float* vert_grid,
	int dem_dim_0, int dem_dim_1, char* geom_type) {

	RTCScene scene = rtcNewScene(device);
  	rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);

  	int num_vert = (dem_dim_0 * dem_dim_1);
  	printf("DEM dimensions: (%d, %d) \n", dem_dim_0, dem_dim_1);
  	printf("Number of vertices: %d \n", num_vert);

	RTCGeometryType rtc_geom_type;
	if (strcmp(geom_type, "triangle") == 0) {
  		rtc_geom_type = RTC_GEOMETRY_TYPE_TRIANGLE;
  	} else if (strcmp(geom_type, "quad") == 0) {
  		rtc_geom_type = RTC_GEOMETRY_TYPE_QUAD;
  	} else { 	
  		rtc_geom_type = RTC_GEOMETRY_TYPE_GRID;
  	}  	

  	RTCGeometry geom = rtcNewGeometry(device, rtc_geom_type);  	
  	rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
  		RTC_FORMAT_FLOAT3, vert_grid, 0, 3*sizeof(float), num_vert);  	
	
	//-------------------------------------------------------------------------
	// Triangle
	//-------------------------------------------------------------------------
	if (strcmp(geom_type, "triangle") == 0) {
		cout << "Selected geometry type: triangle" << endl;
  		int num_tri = ((dem_dim_0 - 1) * (dem_dim_1 - 1)) * 2;
  		printf("Number of triangles: %d \n", num_tri);
  		Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(geom,
  			RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(Triangle),
  			num_tri);
  		int n = 0;
  		for (int i = 0; i < (dem_dim_0 - 1); i++) {
  			for (int j = 0; j < (dem_dim_1 - 1); j++) {
  	  			triangles[n].v0 = (i * dem_dim_1) + j;
  	  			triangles[n].v1 = (i * dem_dim_1) + j + 1;
  	  			triangles[n].v2 = ((i + 1) * dem_dim_1) + j;
	  			n++;
  	  			triangles[n].v0 = (i * dem_dim_1) + j + 1;
  	  			triangles[n].v1 = ((i + 1) * dem_dim_1) + j + 1;
  	  			triangles[n].v2 = ((i + 1) * dem_dim_1) + j;
  	  			n++;
  			}
  		}
	//-------------------------------------------------------------------------
	// Quad
	//-------------------------------------------------------------------------
  	} else if (strcmp(geom_type, "quad") == 0) {
  		cout << "Selected geometry type: quad" << endl;
		int num_quad = ((dem_dim_0 - 1) * (dem_dim_1 - 1));
  		printf("Number of quads: %d \n", num_quad);							   
  		Quad* quads = (Quad*) rtcSetNewGeometryBuffer(geom,
  			RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, sizeof(Quad),
  			num_quad);
  		int n = 0;
  		for (int i = 0; i < (dem_dim_0 - 1); i++) {
  			for (int j = 0; j < (dem_dim_1 - 1); j++) {
  			//  identical to grid scene (-> otherwise reverse v0, v1, ...)
  	  		quads[n].v0 = (i * dem_dim_1) + j;
  	  		quads[n].v1 = (i * dem_dim_1) + j + 1;
  	  		quads[n].v2 = ((i + 1) * dem_dim_1) + j + 1;
  	  		quads[n].v3 = ((i + 1) * dem_dim_1) + j;
  	  		n++;
  		}
  	}    	
	//-------------------------------------------------------------------------
	// Grid
	//-------------------------------------------------------------------------  	
  	} else {
  		cout << "Selected geometry type: grid" << endl;
		RTCGrid* grid = (RTCGrid*)rtcSetNewGeometryBuffer(geom,
			RTC_BUFFER_TYPE_GRID, 0, RTC_FORMAT_GRID, sizeof(RTCGrid), 1);
    	grid[0].startVertexID = 0;
    	grid[0].stride        = dem_dim_1;
    	grid[0].width         = dem_dim_1;
    	grid[0].height        = dem_dim_0;
  	}
	//-------------------------------------------------------------------------

	auto start = std::chrono::high_resolution_clock::now();

	// Commit geometry
	rtcCommitGeometry(geom);

	rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	//-------------------------------------------------------------------------

	// Commit scene
	rtcCommitScene(scene);

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time = end - start;
	cout << "BVH build time: " << time.count() << " s" << endl;

	return scene;

}

//#############################################################################

Rectangle::Rectangle(int X0, int Y0, int X1, int Y1) {
    x0 = X0;
    y0 = Y0;
    x1 = X1;
    y1 = Y1;
    
    device = initializeDevice();
    
}

Rectangle::~Rectangle() {
}

int Rectangle::getArea() {
    return (x1 - x0) * (y1 - y0);
}

void Rectangle::move(int dx, int dy) {
    x0 += dx;
    y0 += dy;
    x1 += dx;
    y1 += dy;
}

void Rectangle::initialise(float* vert_grid,
	int dem_dim_0, int dem_dim_1,
	char* geom_type,
	int offset_0, int offset_1,
	float* vec_tilt,
	float* vec_norm,
	int dim_in_0, int dim_in_1) {

	dem_dim_0_cl = dem_dim_0;
	dem_dim_1_cl = dem_dim_1;
	vert_grid_cl = vert_grid;
	offset_0_cl = offset_0;
	offset_1_cl = offset_1;	
	vec_tilt_cl = vec_tilt;
	vec_norm_cl = vec_norm;
	dim_in_0_cl = dim_in_0;
	dim_in_1_cl = dim_in_1;

	auto start_ini = std::chrono::high_resolution_clock::now();

	scene = initializeScene(device, vert_grid, dem_dim_0, dem_dim_1, geom_type);
	
	auto end_ini = std::chrono::high_resolution_clock::now();
  	std::chrono::duration<double> time = end_ini - start_ini;
  	cout << "Total initialisation time: " << time.count() << " s" << endl;
  	
  	// Allocate array (dim_in_0 * dim_in_1) for surface enlargement factor
  	// surface_enl_fac (constant value; does not chance with moving sun)

}

void Rectangle::shadow(float* sun_position, float* shaddow_buffer) {

	tbb::parallel_for(tbb::blocked_range<size_t>(0,dim_in_0_cl),
		[&](tbb::blocked_range<size_t> r) {  // parallel

	// for (size_t i = 0; i < (size_t)dim_in_0_cl; i++) {  // serial
	for (size_t i=r.begin(); i<r.end(); ++i) {  // parallel
  		for (size_t j = 0; j < (size_t)dim_in_1_cl; j++) {

    		// Get components of terrain surface normal
    		size_t ind_vec = lin_ind_2d(dim_in_1_cl, i, j) * 3;
  			float tilt_x = vec_tilt_cl[ind_vec];
  			ind_vec += 1;
  			float tilt_y = vec_tilt_cl[ind_vec];
  			ind_vec += 1;
  			float tilt_z = vec_tilt_cl[ind_vec];
  
  			// Ray origin
  			size_t ind_2d = lin_ind_2d(dem_dim_1_cl, i + offset_0_cl,
  				j + offset_1_cl);
  			float ray_org_x = vert_grid_cl[ind_2d * 3 + 0];
  			float ray_org_y = vert_grid_cl[ind_2d * 3 + 1];
  			float ray_org_z = vert_grid_cl[ind_2d * 3 + 2] + 1.0; // ---- temporary! not a good solution!!!!

  			// Sun vector
  			float sun_x = (sun_position[0] - ray_org_x);
  			float sun_y = (sun_position[1] - ray_org_y);
  			float sun_z = (sun_position[2] - ray_org_z);
  			float mag = sqrt(sun_x * sun_x + sun_y * sun_y + sun_z * sun_z);
  			sun_x = sun_x / mag;
  			sun_y = sun_y / mag;
  			sun_z = sun_z / mag;
  			
  			// Check for self-shadowing
  			float dot_prod = tilt_x * sun_x + tilt_y * sun_y + tilt_z * sun_z;
  			size_t ind_arr = lin_ind_2d(dim_in_1_cl, i, j);
  			if (dot_prod > 0.0) {
  			
				// Intersect context
  				struct RTCIntersectContext context;
  				rtcInitIntersectContext(&context);

  				// Ray structure
  				struct RTCRay ray;
  				ray.org_x = ray_org_x;
  				ray.org_y = ray_org_y;
  				ray.org_z = ray_org_z;
  				ray.dir_x = sun_x;  // use unit vector later!
  				ray.dir_y = sun_y;
  				ray.dir_z = sun_z;
  				ray.tnear = 0.0;
  				//ray.tfar = std::numeric_limits<float>::infinity();
  				ray.tfar = 100000.0;
  				//ray.mask = -1;
  				//ray.flags = 0;

  				// Intersect ray with scene
  				rtcOccluded1(scene, &context, &ray);
// 				cout << "ray.tfar: " << ray.tfar << endl;
// 				cout << "dem_dim_0: " << dem_dim_0_cl << endl;
// 				cout << "dem_dim_1: " << dem_dim_1_cl << endl;
// 				cout << "vert_grid[8]: " << vert_grid_cl[8] << endl;
// 				cout << "vert_grid[33]: " << vert_grid_cl[33] << endl;


				if (ray.tfar < 0.0) {
					shaddow_buffer[ind_arr] = 2.0;
				} else {
					shaddow_buffer[ind_arr] = 0.0;
				}
			
			} else {
			
				shaddow_buffer[ind_arr] = 1.0;
			
			}
	
		}
	}
	
	}); // parallel

}


