#include "PairwiseRegistration.h"

PairwiseRegistration::PairwiseRegistration()
{}

PairwiseRegistration::PairwiseRegistration(pcl::PointCloud<PointT> &cloud_prev_arg, pcl::PointCloud<PointT> &cloud_curr_arg)
{
	this->cloud_prev = cloud_prev_arg.makeShared();
	this->cloud_curr = cloud_curr_arg.makeShared();
}

PairwiseRegistration::~PairwiseRegistration()
{}

void PairwiseRegistration::updateParameters(pcl::PointCloud<PointT> &cloud_prev_arg, pcl::PointCloud<PointT> &cloud_curr_arg, Eigen::Matrix4f tf_matrix_prev_arg)
{
	this->cloud_prev = cloud_prev_arg.makeShared();
	this->cloud_curr = cloud_curr_arg.makeShared();
	this->tf_matrix_prev = tf_matrix_prev_arg;
}

void PairwiseRegistration::runRegistration()
{
	pcl::IterativeClosestPoint<PointT, PointT> icp;
	icp.setInputSource(this->cloud_curr);
	icp.setInputTarget(this->cloud_prev);
	icp.align(this->cloud_registered);
	this->tf_matrix = icp.getFinalTransformation() * this->tf_matrix_prev;
}

Eigen::Matrix4f PairwiseRegistration::getTfMatrix()
{
	return this->tf_matrix;
}

pcl::PointCloud<PointT> PairwiseRegistration::getRegisteredCloud()
{
	return this->cloud_registered;
}
