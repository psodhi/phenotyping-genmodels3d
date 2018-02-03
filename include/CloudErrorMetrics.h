/**
 * @file CloudErrorMetrics.h
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu)
 * @date 2017-06-15
 * @brief A class to compute a scalar hausdorff error value between source and target point clouds
 */

#ifndef CloudErrorMetrics_H_
#define CloudErrorMetrics_H_

#include "pcl-utils.h"

using namespace std;

class CloudErrorMetrics {

public:

	CloudErrorMetrics();
	CloudErrorMetrics(pcl::PointCloud<PointT> &, pcl::PointCloud<PointT> &);
	~CloudErrorMetrics();

	void updateParameters(pcl::PointCloud<PointT> &, pcl::PointCloud<PointT> &);
	float getErrorValue();
	float compute_hausdorff_distance();

private:
	pcl::PointCloud<PointT>::Ptr cloud_src, cloud_tgt;
	float error_val;
};

#endif