/**
 * @file PhytomerGeneration.h
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu)
 * @date 2017-06-01
 * @brief A class to generate a parametric 3D plant phytomer using vtk/pcl libraries
 */

#ifndef PhytomerGeneration_H_
#define PhytomerGeneration_H_

#include "PointsFromSurfaceSampler.h"

using namespace std;

// endverts hold start and end vertices 
struct endverts {
    Eigen::Vector3f startpos;
    Eigen::Vector3f endpos;
};

struct pose {
	Eigen::Vector3f position;
	Eigen::Vector3f orientation;
};


class PhytomerGeneration {

public:

	PhytomerGeneration();
	~PhytomerGeneration();

	void updateParameters(float stem_radius_arg, float leaf_angle_arg, float leaf_length_arg);
	void updateParameters(float stem_radius_arg, float leaf_angle_arg, float leaf_length_arg, float stem_internode_length_arg, float leaf_rotate_cont_angle_arg, float leaf_width_arg);
	void setNumberOfSamplePoints(int sample_points_arg);
	pcl::PointCloud<PointT> getPhytomerPointCloudNoisy(float noise_perc);
	Eigen::Matrix3f make_rotation_matrix(float leaf_angle_rad, char symbol);
	float deg2rad(float);
	float fluctuate(float, float);
	tuple<Eigen::Matrix<float, 4, 3> , Eigen::Vector3f> get_leaf_cross_section(Eigen::Vector3f, int, Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f, int);
	void make_leaf_patches(vector<endverts>, int, vtkSmartPointer<vtkPolyData> &);
	void plot_leaves(vector<endverts>, vtkSmartPointer<vtkPolyData> &, pcl::PointCloud<PointT>::Ptr);
	void plot_stem_nodes(vector<endverts>, vtkSmartPointer<vtkPolyData> &, pcl::PointCloud<PointT>::Ptr);
	void plot_from_string();
	void generate_sorghum_phytomer();
	pcl::PointCloud<PointT> getPhytomerPointCloud();

	void colorizeCloud(pcl::PointCloud<PointT> &, int[]);
	int visualizeVTKPolyData(vtkSmartPointer<vtkPolyData>);
	int visualizeVTKPolyData(vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPolyData>);

	int main();

private:

	string lsys_string_;
	int sample_points_;
	float percentage_stem_points_;
	int scale_;

	float stem_radius_;
	float stem_internode_length_;

	float leaf_angle_;
	float leaf_length_;
	float leaf_width_;
	float leaf_thickness_;
	float leaf_rotate_cont_angle_;
	float leaf_area_;

	pcl::PointCloud<PointT> cloud_gen_phytomer_;

};
  
#endif