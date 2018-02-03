#include "CrossEntropyOptimization.h"

using namespace std;

CrossEntropyOptimization::CrossEntropyOptimization() {
	vis_flag_ = true;
}

CrossEntropyOptimization::CrossEntropyOptimization(std::string cfg_filename) {
	cfg_filename_ = cfg_filename;
	vis_flag_ = true;
	readParamsFromFile();
}

void CrossEntropyOptimization::setVisualizationFlag(bool vis_flag) {
	vis_flag_ = vis_flag;
}

void CrossEntropyOptimization::readParamsFromFile() {

	// Read config parameters from cfg file passed as argument

	std::vector<std::string> cfg_params;
	std::ifstream cfg_ifstream(cfg_filename_);
	 
	if(cfg_ifstream.is_open()) {
	   std::string cfg_line;

	   std::cout << "Input Set of Parameters Read : " << std::endl;
	   while(std::getline(cfg_ifstream, cfg_line)) {
	       if(cfg_line[0] == '$') {
	       	std::cout << cfg_line << " : ";
	       	std::getline(cfg_ifstream, cfg_line);
	       	std::cout << cfg_line << std::endl;
			cfg_params.push_back(cfg_line);
	      }
	  }
	}
	else {
		std::cout << "ERROR opening config parameters file" << std::endl;
	}

	// Initialize parameters using params read from the cfg file

	dir_src_phytomers_ = cfg_params[0];
	dir_dest_ceopt_ = cfg_params[1];

	num_samples_ = std::stoi(cfg_params[2]); // number of samples (or particles) in each iteration
	percentile_elite_ = std::stod(cfg_params[3]); // elite sample percentile
	num_opt_iters_ = std::stoi(cfg_params[4]); // number of optimization iterations
	noise_factor_sigma_ = std::stod(cfg_params[5]); 

	mu_stemdiameter_ = std::stod(cfg_params[6]);
	mu_leafangle_ = std::stod(cfg_params[7]);
	mu_leaflength_ = std::stod(cfg_params[8]);
	mu_steminternodelen_ = std::stod(cfg_params[9]); // mu_leaflength/6 + 1.5
	mu_leafrotateang_ = std::stod(cfg_params[10]);
	mu_leafwidth_ = std::stod(cfg_params[11]); // mu_leaflength/10 + 1.5

	sigma_stemdiameter_ = std::stod(cfg_params[12]);
	sigma_leafangle_ = std::stod(cfg_params[13]);
	sigma_leaflength_ = std::stod(cfg_params[14]);
	sigma_steminternodelen_ = std::stod(cfg_params[15]);
	sigma_leafrotateang_ = std::stod(cfg_params[16]);
	sigma_leafwidth_ = std::stod(cfg_params[17]);

	// Compute some additional optimization related params
	num_optvars_ = 6; // number of unknown phenotype variables being optimized
	num_samples_elite_ = ceil((1-percentile_elite_)*num_samples_); // use ceil over floor to avoid num_samples_elite becoming 0
}

// function that returns sorted indices
template <typename T>
inline std::vector<size_t> ordered(std::vector<T> const& values)
{
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

template <typename T>
inline vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  vector<T> xs(N);
  typename vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}

void CrossEntropyOptimization::start()
{
	// Set point cloud related variables 
	pcl::PointCloud<PointT>::Ptr cloud_phytomer_ref (new pcl::PointCloud<PointT>);
	pcl::PointCloud<PointT>::Ptr cloud_phytomer_sample (new pcl::PointCloud<PointT>);
	Eigen::Vector4f pt_centroid;

	// Set ICP & Hausdorff related variables
	Eigen::Matrix4f tf_matrix, tf_matrix_prev;
	tf_matrix_prev << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;
	float err_hausdorff;

	// Set dir/files related variables 
	std::string filename;
	std::stringstream ss_dir;

	vector<::boost::filesystem::path> files_src;
	std::string delimiter = ".";
	bool isPLY = false;
	pclutils::get_ext_dir(dir_src_phytomers_, ".pcd", files_src);

	for (int i = 0; i < files_src.size(); i++) {

		// Reset parameters for each new reference phytomer point cloud whose phenotypes need to be estimated
		readParamsFromFile();

		// Read reference phytomer point cloud whose phenotypes need to be estimated
		filename = files_src[i].string();
		filename = filename.substr(0, filename.find(delimiter));

		cout << "Reference Phytomer File Name : " << filename << endl;
		ss_dir.str(""); ss_dir << dir_src_phytomers_.string() << files_src[i].string();

		if (isPLY == true) {
			if (pcl::io::loadPLYFile(ss_dir.str(), *cloud_phytomer_ref) < 0 ) {
				cerr << "Cloud Reading Failed. Skipping this Cloud." << endl;
				continue;
			}
		}
		else {
			if (pcl::io::loadPCDFile(ss_dir.str(), *cloud_phytomer_ref) < 0 ) {
				cerr << "Cloud Reading Failed. Skipping this Cloud." << endl; 
				continue;
			}
		}

		pclutils::filterCloudNaNs(cloud_phytomer_ref, cloud_phytomer_ref);

		Eigen::VectorXf mean_vec(6);
		mean_vec(0) = mu_stemdiameter_; mean_vec(1) = mu_leafangle_, mean_vec(2) = mu_leaflength_;
		mean_vec(3) = mu_steminternodelen_, mean_vec(4) = mu_leafrotateang_, mean_vec(5) = mu_leafwidth_;

		Eigen::MatrixXf covariance(num_optvars_, num_optvars_);
		covariance << pow(sigma_stemdiameter_,1), 0, 0, 0, 0, 0,
						0, pow(sigma_leafangle_,1), 0, 0, 0, 0, 
						0, 0, pow(sigma_leaflength_,1), 0, 0, 0,
						0, 0, 0, pow(sigma_steminternodelen_,1), 0, 0,
						0, 0, 0, 0, pow(sigma_leafrotateang_,1), 0,
						0, 0, 0, 0, 0, pow(sigma_leafwidth_,1);

		vector<float> noise_std_stemdiameter = linspace(noise_factor_sigma_*sigma_stemdiameter_, 0.0f, num_opt_iters_);
		vector<float> noise_std_leafangle = linspace(noise_factor_sigma_*sigma_leafangle_, 0.0f, num_opt_iters_);
		vector<float> noise_std_leaflength = linspace(noise_factor_sigma_*sigma_leaflength_, 0.0f, num_opt_iters_);
		vector<float> noise_std_steminternodelen = linspace(noise_factor_sigma_*sigma_steminternodelen_, 0.0f, num_opt_iters_);
		vector<float> noise_std_leafrotateang = linspace(noise_factor_sigma_*sigma_leafrotateang_, 0.0f, num_opt_iters_);
		vector<float> noise_std_leafwidth = linspace(noise_factor_sigma_*sigma_leafwidth_, 0.0f, num_opt_iters_);

		cout << " ====== Initial Values " << " ======" << endl;
		cout << "Mean (stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width) : " << endl << mean_vec << endl;
		cout << "Covariance (stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width) : " << endl << covariance << endl;

		pcl::visualization::PCLVisualizer viewer1("cloud_phytomer");
		viewer1.setBackgroundColor (255.0, 255.0, 255.0);

		for (int iter = 0; iter < num_opt_iters_; iter++) {

			// Visualize phytomer generated using mean value of all samples aligned against the reference phytomer
			if (vis_flag_) {
				viewer1.removeAllPointClouds();
				// viewer1.loadCameraParameters("top_side_view_phytomer_greenhouse.cam");

				PairwiseRegistration reg1;
				PhytomerGeneration phygen1;
				std::stringstream ss_name;

				pcl::PointCloud<PointT> cloud_phytomer_vis;
				phygen1.updateParameters(0.5*mean_vec(0), mean_vec(1), mean_vec(2), mean_vec(3), mean_vec(4), mean_vec(5));
				phygen1.generate_sorghum_phytomer();
				cloud_phytomer_vis = phygen1.getPhytomerPointCloud();

				compute3DCentroid (cloud_phytomer_vis, pt_centroid);
				demeanPointCloud (cloud_phytomer_vis, pt_centroid, cloud_phytomer_vis);

				reg1.updateParameters(*cloud_phytomer_ref, cloud_phytomer_vis, tf_matrix_prev);
				reg1.runRegistration();
				pcl::transformPointCloud(cloud_phytomer_vis, cloud_phytomer_vis, reg1.getTfMatrix());

				pcl::visualization::PointCloudColorHandlerCustom<PointT> colorHandlerPhytomerRef(cloud_phytomer_ref, 0, 0, 0);
				ss_name.str(""); ss_name << "cloud_ref_phytomer_" << iter;
				viewer1.addPointCloud(cloud_phytomer_ref, colorHandlerPhytomerRef, ss_name.str());
				viewer1.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1.5, ss_name.str());
				ss_name.str(""); ss_name << "cloud_genmean_phytomer_" << iter;
				viewer1.addPointCloud(cloud_phytomer_vis.makeShared(), ss_name.str());
				viewer1.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1.5, ss_name.str());

				// viewer1.addCoordinateSystem();
				// viewer1.setSize(1920, 1080);
				viewer1.spinOnce();

				ss_name.str(""); ss_name << dir_dest_ceopt_.string() << "/" << filename << "_iter_" << setfill('0') << setw(2) << iter << ".png";
				viewer1.saveScreenshot(ss_name.str());
			}
			
			// Cross entropy importance sampling

			Eigen::EigenMultivariateNormal<float> mvnObj(mean_vec, covariance);
			auto samples_vals = mvnObj.samples(num_samples_); // sampled_vals is (num_optvars x num_samples_) sized Eigen::Matrix
			vector<float> error_vector(num_samples_); 

			vector<pcl::PointCloud<PointT> > cloud_phytomer_vector(num_samples_);
			omp_set_num_threads(8);
			#pragma omp parallel for
			for (int i = 0; i < num_samples_; i++)
			{
				// cout << "Particle ID : " << i << endl;
				PairwiseRegistration reg1;
				CloudErrorMetrics cloudErr1;
				PhytomerGeneration phygen1;
				Eigen::Vector4f pt_centroid_particle;

				// Phytomer Generation : updateParameters(stem_radius, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width)
				phygen1.setNumberOfSamplePoints(cloud_phytomer_ref->points.size());
				phygen1.updateParameters( std::max( abs(0.5f*samples_vals(0, i)), 1.5f), std::max( abs(samples_vals(1, i)), 10.0f), std::max( abs(samples_vals(2, i)), 5.0f),
											std::max( abs(samples_vals(3, i)), 2.0f), std::max( abs(samples_vals(4, i)), 1.5f), std::max( abs(samples_vals(5, i)), 1.0f) );
				phygen1.generate_sorghum_phytomer();
				cloud_phytomer_vector[i] = phygen1.getPhytomerPointCloud();

				compute3DCentroid (cloud_phytomer_vector[i], pt_centroid_particle);
				demeanPointCloud (cloud_phytomer_vector[i], pt_centroid_particle, cloud_phytomer_vector[i]);

				// ICP Registration : updateParameters(cloud_tgt, cloud_src, tf_matrix_prev)
				reg1.updateParameters(*cloud_phytomer_ref, cloud_phytomer_vector[i], tf_matrix_prev);
				reg1.runRegistration();
				pcl::transformPointCloud(cloud_phytomer_vector[i], cloud_phytomer_vector[i], reg1.getTfMatrix());

				// // Hausdorff Distance Error : updateParameters(cloud_ref, cloud_gen) 
				cloudErr1.updateParameters(*cloud_phytomer_ref, cloud_phytomer_vector[i]);
				error_vector[i] = cloudErr1.compute_hausdorff_distance();
			}

			// Update Mean/Covariance using elite set of samples
			// sort score_vector and get sort_idx and choose top num_samples_elite_ in sort_idx & recompute mean, sigma.
			std::vector<size_t> error_sort_idxs = ordered(error_vector);

			// accumulate elite sample values & compute mean values
			mean_vec(0) = 0.0f; mean_vec(1) = 0.0f; mean_vec(2) = 0.0f; mean_vec(3) = 0.0f; mean_vec(4) = 0.0f; mean_vec(5) = 0.0f;
			for (int i = 0; i < num_samples_elite_; i++)
			{
				mean_vec(0) = mean_vec(0) + std::abs(samples_vals(0, error_sort_idxs[i]));
				mean_vec(1) = mean_vec(1) + std::abs(samples_vals(1, error_sort_idxs[i]));
				mean_vec(2) = mean_vec(2) + std::abs(samples_vals(2, error_sort_idxs[i]));
				mean_vec(3) = mean_vec(3) + std::abs(samples_vals(3, error_sort_idxs[i]));
				mean_vec(4) = mean_vec(4) + std::abs(samples_vals(4, error_sort_idxs[i]));
				mean_vec(5) = mean_vec(5) + std::abs(samples_vals(5, error_sort_idxs[i]));
			}
			mean_vec(0) = mean_vec(0)/num_samples_elite_; mean_vec(1) = mean_vec(1)/num_samples_elite_; mean_vec(2) = mean_vec(2)/num_samples_elite_;
			mean_vec(3) = mean_vec(3)/num_samples_elite_; mean_vec(4) = mean_vec(4)/num_samples_elite_; mean_vec(5) = mean_vec(5)/num_samples_elite_;

			// update covariance matrix values using elite sample values
			covariance(0,0) = 0.0f; covariance(1,1) = 0.0f; covariance(2,2) = 0.0f;
			covariance(3,3) = 0.0f; covariance(4,4) = 0.0f; covariance(5,5) = 0.0f;
			for (int i = 0; i < num_samples_elite_; i++)
			{
				covariance(0,0) = covariance(0,0) + std::pow( std::abs(samples_vals(0, error_sort_idxs[i]))-mean_vec(0), 2);
				covariance(1,1) = covariance(1,1) + std::pow( std::abs(samples_vals(1, error_sort_idxs[i]))-mean_vec(1), 2);
				covariance(2,2) = covariance(2,2) + std::pow( std::abs(samples_vals(2, error_sort_idxs[i]))-mean_vec(2), 2);
				covariance(3,3) = covariance(3,3) + std::pow( std::abs(samples_vals(3, error_sort_idxs[i]))-mean_vec(3), 2);
				covariance(4,4) = covariance(4,4) + std::pow( std::abs(samples_vals(4, error_sort_idxs[i]))-mean_vec(4), 2);
				covariance(5,5) = covariance(5,5) + std::pow( std::abs(samples_vals(5, error_sort_idxs[i]))-mean_vec(5), 2);
			}

			covariance(0,0) = std::sqrt(covariance(0,0))/num_samples_elite_ + std::pow(noise_std_stemdiameter[iter], 1);
			covariance(1,1) = std::sqrt(covariance(1,1))/num_samples_elite_ + std::pow(noise_std_leafangle[iter], 1);
			covariance(2,2) = std::sqrt(covariance(2,2))/num_samples_elite_ + std::pow(noise_std_leaflength[iter], 1);
			covariance(3,3) = std::sqrt(covariance(3,3))/num_samples_elite_ + std::pow(noise_std_steminternodelen[iter], 1);
			covariance(4,4) = std::sqrt(covariance(4,4))/num_samples_elite_ + std::pow(noise_std_leafrotateang[iter], 1);
			covariance(5,5) = std::sqrt(covariance(5,5))/num_samples_elite_ + std::pow(noise_std_leafwidth[iter], 1);

			cout << " ====== Optimization Iteration : " << iter+1 << " ======" << endl;
			cout << "Mean (stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width) : " << endl << mean_vec << endl;
			cout << "Covariance (stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width) : " << endl << covariance << endl;

			// Write mean/stdev output to file at end of each iteration 
			ss_dir.str(""); ss_dir << dir_dest_ceopt_.string() << "/" << filename << ".txt";
			cout << "Writing output to file : " << ss_dir.str() << endl;
			ofstream outfile;
			outfile.open(ss_dir.str());
			outfile << "Order of Variables : stem_diameter, leaf_angle, leaf_length, stem_internode_len, leaf_rotate_angle, leaf_width \n"; 
			outfile << "Mean & Standard Deviations after " <<  (iter+1) << " iterations for " << filename << "\n\n";
			outfile << mean_vec << "\n\n";
			outfile << covariance << "\n";
			outfile.close();
		}

		cloud_phytomer_ref->clear();
	}
}