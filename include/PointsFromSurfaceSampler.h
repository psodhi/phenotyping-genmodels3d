/**
 * @file PointsFromSurfaceSampler.h
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu) (functions adapted from pcl/tools/mesh_sampling.cpp)
 * @date 2017-06-10
 * @brief A class to uniformly sample 3D points from a surface mesh
 */

#ifndef PointsFromSurfaceSampler_H_
#define PointsFromSurfaceSampler_H_

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <chrono>

// PCL related include files
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/io/vtk_lib_io.h>

#include <pcl/impl/point_types.hpp>
#include <pcl/PolygonMesh.h>

#include <pcl/common/common.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>

#include <pcl/visualization/pcl_visualizer.h>

// VTK related include files
#include <vtkCylinderSource.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

using namespace std;
typedef pcl::PointXYZRGB PointT;

class PointsFromSurfaceSampler {

public:

	PointsFromSurfaceSampler();
	PointsFromSurfaceSampler(vtkSmartPointer<vtkPolyData>);
	~PointsFromSurfaceSampler();

	double uniform_deviate (int);
	void randomPointTriangle (float, float, float, float, float, float, float, float, float, Eigen::Vector4f&);
	void randPSurface (vtkPolyData*, std::vector<double>* , double, Eigen::Vector4f&, bool, Eigen::Vector3f&);
	void uniform_sampling (vtkSmartPointer<vtkPolyData>, size_t, bool, pcl::PointCloud<pcl::PointNormal>&);
	void uniform_sampling (vtkSmartPointer<vtkPolyData>, size_t, bool, pcl::PointCloud<PointT>&);

private:
	vtkSmartPointer<vtkPolyData> polydata_;
};
  
#endif