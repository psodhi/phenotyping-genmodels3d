#include "CrossEntropyOptimization.h"
#include "PhytomerGeneration.h"
#include "cmaes.h"

using namespace std;
using namespace libcmaes;

/**
 * @file main.cpp
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu)
 * @date 2017-06-01
 * @brief main.cpp calling functions to generate phytomers or call cross-entropy optimization
 */

void generatePhytomer()
{
	float stem_diameter = 3; // cm
	float leaf_angle = 45; // degrees
	float leaf_length = 50; // cm
	float stem_internode_len = 10;
	float leaf_rotate_angle = 3;
	float leaf_width = 5;

	pcl::PointCloud<PointT> cloud_gen_phytomer;

	PhytomerGeneration phygen1;
	phygen1.updateParameters(0.5f*stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width);
	phygen1.setNumberOfSamplePoints(10000);
	phygen1.plot_from_string();
	
	cloud_gen_phytomer = phygen1.getPhytomerPointCloud();
	// cloud_gen_phytomer = phygen1.getPhytomerPointCloudNoisy(0.05f); // 5% gaussian noise

	// Visualization of generated plant phytomer point cloud
	pcl::visualization::PCLVisualizer viewer1("cloud_gen_phytomer");
	viewer1.addPointCloud(cloud_gen_phytomer.makeShared());
	viewer1.setBackgroundColor (1.0, 1.0, 1.0);
	viewer1.addCoordinateSystem ();
	while (!(viewer1.wasStopped()))
		viewer1.spinOnce();
}


int main(int argc, char** argv)
{
	std::string cfg_filename = argv[1];
	CrossEntropyOptimization ceopt_obj(cfg_filename);
	ceopt_obj.start();

	return 0;
}