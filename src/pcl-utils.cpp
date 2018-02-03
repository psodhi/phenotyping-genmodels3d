/**
* @file pcl-utils.cpp
* @author Paloma Sodhi (psodhi@andrew.cmu.edu)
* @date 2017-01-15
* @brief A set of utility functions to work with the PCL library.
* (implementation file)
*/

#include "pcl-utils.h"

namespace pclutils { 

	std::vector <std::string> read_directory( const std::string& path = std::string() )
	{
	  std::vector <std::string> result;
	  dirent* de;
	  DIR* dp;
	  errno = 0;
	  dp = opendir( path.empty() ? "." : path.c_str() );
	  if (dp)
	    {
	    while (true)
	      {
	      errno = 0;
	      de = readdir( dp );
	      if (de == NULL) break;
	      result.push_back( std::string( de->d_name ) );
	      }
	    closedir( dp );
	    std::sort( result.begin(), result.end() );
	    }
	  return result;
	}

	void get_ext_dir(const ::boost::filesystem::path& root, const string& ext, vector<::boost::filesystem::path>& ret)
	{
	    if(!::boost::filesystem::exists(root) || !::boost::filesystem::is_directory(root)) return;
	    ::boost::filesystem::recursive_directory_iterator it(root);
	    ::boost::filesystem::recursive_directory_iterator endit;
	    while(it != endit)
	    {
	        if(::boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) ret.push_back(it->path().filename());
	        ++it;
	    }
	    std::sort(ret.begin(), ret.end());
	}

	void filterCloudNaNs(pcl::PointCloud<PointT>::Ptr cloud_in, pcl::PointCloud<PointT>::Ptr cloud_out)
	{
		  std::vector<int> indices;
		  pcl::removeNaNFromPointCloud(*cloud_in, *cloud_out, indices);
	}


	void voxelizeCloud(pcl::PointCloud<PointT>::Ptr cloud_in, float leaf_size, pcl::PointCloud<PointT>::Ptr cloud_out)
	{
		pcl::VoxelGrid<PointT> vor;
		vor.setInputCloud (cloud_in);
		vor.setLeafSize (leaf_size, leaf_size, leaf_size);
		vor.filter (*cloud_out);
	}

	void passThroughFilter(pcl::PointCloud<PointT>::Ptr cloud_in, float min[], float max[], pcl::PointCloud<PointT>::Ptr cloud_out)
	{
		pcl::PassThrough<PointT> pass;
		
		pass.setInputCloud(cloud_in);
		pass.setFilterFieldName("x");
		pass.setFilterLimits(min[0], max[0]);
		pass.filter(*cloud_out);

		pass.setInputCloud(cloud_out);
		pass.setFilterFieldName("y");
		pass.setFilterLimits(min[1], max[1]);
		pass.filter(*cloud_out);
		
		pass.setInputCloud(cloud_out);
		pass.setFilterFieldName("z");
		pass.setFilterLimits(min[2], max[2]);
		pass.filter(*cloud_out);

	}

	void statisticalOutlierFilter(pcl::PointCloud<PointT>::Ptr cloud_in, float mean, float std_dev, pcl::PointCloud<PointT>::Ptr cloud_out)
	{
		pcl::StatisticalOutlierRemoval<PointT> sor;
		
		sor.setInputCloud(cloud_in);
		sor.setMeanK(mean);
		sor.setStddevMulThresh(std_dev);
		sor.filter(*cloud_out);
	}

	void transformCloudXYZ(pcl::PointCloud<PointT>::Ptr cloud_in, float trans[], float rot[], pcl::PointCloud<PointT>::Ptr cloud_out)
	{
		Eigen::Affine3f transform = Eigen::Affine3f::Identity();
		transform.translation() << trans[0], trans[1], trans[2];
		transform.rotate (Eigen::AngleAxisf (rot[0], Eigen::Vector3f::UnitX()));
		transform.rotate (Eigen::AngleAxisf (rot[1], Eigen::Vector3f::UnitY()));
		transform.rotate (Eigen::AngleAxisf (rot[2], Eigen::Vector3f::UnitZ()));

		pcl::transformPointCloud (*cloud_in, *cloud_out, transform);
	}

	void colorizeCloud(pcl::PointCloud<PointT>::Ptr cloud_in, int rgb_vec[])
	{
		for (int i = 0; i < cloud_in->points.size(); i++)
		{
			cloud_in->points[i].r = rgb_vec[0];
			cloud_in->points[i].g = rgb_vec[1];
			cloud_in->points[i].b = rgb_vec[2];
		}
	}

	void visualizeSingleCloud(pcl::PointCloud<PointT>::Ptr cloud_in, string window_name)
	{
		pcl::visualization::PCLVisualizer viewer1(window_name);
		viewer1.addPointCloud(cloud_in, window_name);
		// viewer1.setBackgroundColor (153/255.0, 204/255.0, 255/255.0);
		viewer1.setBackgroundColor (1.0, 1.0, 1.0);
		viewer1.addCoordinateSystem ();
		while (!(viewer1.wasStopped()))
			viewer1.spinOnce();
	}

	void visualizeTwoCloudsSingleViewer(pcl::PointCloud<PointT>::Ptr cloud_in1, pcl::PointCloud<PointT>::Ptr cloud_in2, string win_name)
	{
		stringstream ss_winname;

		pcl::visualization::PCLVisualizer viewer1(win_name);
		viewer1.setBackgroundColor (153/255.0, 204/255.0, 255/255.0);

		ss_winname.str(""); ss_winname << win_name << "_1";
		viewer1.addPointCloud(cloud_in1, ss_winname.str());
		ss_winname.str(""); ss_winname << win_name << "_2";
		viewer1.addPointCloud(cloud_in2, ss_winname.str());
		viewer1.addCoordinateSystem();

		while ( !(viewer1.wasStopped()) ) {
			viewer1.spin();
		}
	}

	void visualizeTwoCloudsTwoViewers(pcl::PointCloud<PointT>::Ptr cloud_in1, pcl::PointCloud<PointT>::Ptr cloud_in2, string win_name1, string win_name2)
	{
		pcl::visualization::PCLVisualizer viewer1(win_name1), viewer2(win_name2);
		pcl::visualization::PointCloudColorHandlerCustom<PointT> colorHandler2(cloud_in2, 255, 0, 0);

		viewer1.setBackgroundColor (153/255.0, 204/255.0, 255/255.0);
		viewer2.setBackgroundColor (153/255.0, 204/255.0, 255/255.0);

		viewer1.addPointCloud(cloud_in1, win_name1);
		viewer2.addPointCloud(cloud_in2, colorHandler2, win_name2);

		while ( !(viewer1.wasStopped()) || !(viewer2.wasStopped()) ) {
			viewer1.spinOnce();
			viewer2.spinOnce();
		}
	}

}