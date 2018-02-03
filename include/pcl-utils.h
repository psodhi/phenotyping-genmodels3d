/**
* @file pcl-utils.h
* @author Paloma Sodhi (psodhi@andrew.cmu.edu)
* @date 2017-01-15
* @brief A set of utility functions to work with the PCL library.
* (header file)
*/

#ifndef PCLUTILS_H_
#define PCLUTILS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>

#include <boost/filesystem.hpp>

#include <pcl/point_types.h>

#include <pcl/common/common.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>

#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/statistical_outlier_removal.h>

#include <pcl/visualization/pcl_visualizer.h>

using namespace std;
typedef pcl::PointXYZRGB PointT;

namespace pclutils {

	std::vector <std::string> read_directory(const std::string&);
	void get_ext_dir(const ::boost::filesystem::path&, const string&, vector<::boost::filesystem::path>&);

	void filterCloudNaNs(pcl::PointCloud<PointT>::Ptr, pcl::PointCloud<PointT>::Ptr);
	void voxelizeCloud(pcl::PointCloud<PointT>::Ptr, float, pcl::PointCloud<PointT>::Ptr);
	void passThroughFilter(pcl::PointCloud<PointT>::Ptr, float[], float[], pcl::PointCloud<PointT>::Ptr);
	void statisticalOutlierFilter(pcl::PointCloud<PointT>::Ptr, float, float, pcl::PointCloud<PointT>::Ptr);
	void transformCloudXYZ(pcl::PointCloud<PointT>::Ptr, float[], float[], pcl::PointCloud<PointT>::Ptr);
	void colorizeCloud(pcl::PointCloud<PointT>::Ptr, int[]);


	void visualizeSingleCloud(pcl::PointCloud<PointT>::Ptr, string);
	void visualizeTwoCloudsSingleViewer(pcl::PointCloud<PointT>::Ptr, pcl::PointCloud<PointT>::Ptr, string);
	void visualizeTwoCloudsTwoViewers(pcl::PointCloud<PointT>::Ptr, pcl::PointCloud<PointT>::Ptr, string, string);
}

#endif