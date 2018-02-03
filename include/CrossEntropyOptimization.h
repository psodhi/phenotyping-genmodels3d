/**
 * @file CrossEntropyOptimization.h
 * @author Paloma Sodhi (psodhi@andrew.cmu.edu)
 * @date 2017-06-20
 * @brief A class employing cross-entropy optimization to resample phytomer model parameter distributions to match reference phytomer 
 */

#ifndef CrossEntropyOptimization_H_
#define CrossEntropyOptimization_H_

#include "pcl-utils.h"
#include "PairwiseRegistration.h"
#include "EigenMultivariateNormal.h"
#include "CloudErrorMetrics.h"
#include "PhytomerGeneration.h"

#include <omp.h>

using namespace std;

class CrossEntropyOptimization {

public:

	CrossEntropyOptimization();
	CrossEntropyOptimization(std::string);
	~CrossEntropyOptimization() {};

	void setVisualizationFlag(bool);
	void readParamsFromFile();
	void start();

private:

	std::string cfg_filename_;
	bool vis_flag_;

	::boost::filesystem::path dir_src_phytomers_;
	::boost::filesystem::path dir_dest_ceopt_;

	int num_samples_; // number of samples (or particles) in each iteration
	float percentile_elite_; // elite sample percentile
	int num_opt_iters_; // number of optimization iterations
	float noise_factor_sigma_; 

	float mu_stemdiameter_;
	float mu_leafangle_;
	float mu_leaflength_;
	float mu_steminternodelen_;
	float mu_leafrotateang_;
	float mu_leafwidth_;

	float sigma_stemdiameter_;
	float sigma_leafangle_;
	float sigma_leaflength_;
	float sigma_steminternodelen_;
	float sigma_leafrotateang_;
	float sigma_leafwidth_;

	int num_optvars_;
	int num_samples_elite_;
};

#endif