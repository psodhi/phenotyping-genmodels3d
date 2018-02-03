/**
 * @file PairwiseRegistration.h
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu)
 * @date 2017-02-01
 * @brief A class to perform 3D registration (ICP) for a point cloud pair
 */

#ifndef PairwiseRegistration_H_
#define PairwiseRegistration_H_

#include "pcl-utils.h"
#include <pcl/registration/icp.h>

using namespace std;
typedef pcl::PointXYZRGB PointT;

class PairwiseRegistration {

public:
	PairwiseRegistration();
	PairwiseRegistration(pcl::PointCloud<PointT> &, pcl::PointCloud<PointT> &);
	~PairwiseRegistration();

	void updateParameters(pcl::PointCloud<PointT> &, pcl::PointCloud<PointT> &, Eigen::Matrix4f);
	void runRegistration();
	Eigen::Matrix4f getTfMatrix();
	pcl::PointCloud<PointT> getRegisteredCloud();

private:
	pcl::PointCloud<PointT>::Ptr cloud_prev, cloud_curr;
	Eigen::Matrix4f tf_matrix_prev, tf_matrix;
	pcl::PointCloud<PointT> cloud_registered;
};

#endif