#include "CloudErrorMetrics.h"

using namespace std;

CloudErrorMetrics::CloudErrorMetrics() {
  this->error_val = 0.0;
}

CloudErrorMetrics::CloudErrorMetrics(pcl::PointCloud<PointT> &cloud_src_arg, pcl::PointCloud<PointT> &cloud_tgt_arg) {
  this->cloud_src = cloud_src_arg.makeShared();
  this->cloud_tgt = cloud_tgt_arg.makeShared();
  this->error_val = 0.0;
}

CloudErrorMetrics::~CloudErrorMetrics()
{}

void CloudErrorMetrics::updateParameters(pcl::PointCloud<PointT> &cloud_src_arg, pcl::PointCloud<PointT> &cloud_tgt_arg) {
  this->cloud_src = cloud_src_arg.makeShared();
  this->cloud_tgt = cloud_tgt_arg.makeShared();
}

float CloudErrorMetrics::getErrorValue() {
  return this->error_val;
}

float CloudErrorMetrics::compute_hausdorff_distance() {

  // compare A:cloud_src to B:cloud_tgt
  pcl::search::KdTree<PointT> tree_b;
  tree_b.setInputCloud (this->cloud_tgt);
  float max_dist_a = -std::numeric_limits<float>::max ();
  for (size_t i = 0; i < (this->cloud_src)->points.size (); ++i)
  {
    std::vector<int> indices (1);
    std::vector<float> sqr_distances (1);

    tree_b.nearestKSearch ((this->cloud_src)->points[i], 1, indices, sqr_distances);
    if (sqr_distances[0] > max_dist_a)
      max_dist_a = sqr_distances[0];
  }

  // compare B:cloud_tgt to A:cloud_src
  pcl::search::KdTree<PointT> tree_a;
  tree_a.setInputCloud (this->cloud_src);
  float max_dist_b = -std::numeric_limits<float>::max ();
  for (size_t i = 0; i < (this->cloud_tgt)->points.size (); ++i)
  {
    std::vector<int> indices (1);
    std::vector<float> sqr_distances (1);

    tree_a.nearestKSearch ((this->cloud_tgt)->points[i], 1, indices, sqr_distances);
    if (sqr_distances[0] > max_dist_b)
      max_dist_b = sqr_distances[0];
  }

  max_dist_a = std::sqrt (max_dist_a);
  max_dist_b = std::sqrt (max_dist_b);

  float dist = std::max (max_dist_a, max_dist_b);

  // DEBUG : Print out error values
  // cout << "A -> B : " << max_dist_a;
  // cout << "   B -> A : " << max_dist_b;
  // cout << "   Hausdorff Distance : " << dist;
  // cout << endl;

  return dist;

}