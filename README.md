# README #

**Generative 3D Models for Estimating Plant Physical Characteristics**

### What is this repository for? ###

This repository contains code for generative model-based phenotyping from sorghum plant phytomers. The approach involves sampling parameterized 3D plant phytomer models from an underlying probability distribution and comparing those against the reference model whose phenotypes we wish to estimate. It then optimizes for an objective of making the mass of this probability distribution approach the true parameters of the reference model.

More details can be found under model-based methods of the following document, 

[1] In-field Plant Phenotyping using Model-free and Model-based methods

<https://www.ri.cmu.edu/publications/in-field-plant-phenotyping-using-model-free-and-model-based-methods/>

### How do I get set up? ###

This code has been tested on a Linux (Ubuntu 16.04) system. 

**Dependencies**

pcl 1.7

Download and install instructions:

<http://pointclouds.org/downloads/>

libvtk6

Download and install instructions:

<https://www.vtk.org/download/>

OpenMP

Check if your compiler supports and implements OpenMP here,

<http://www.openmp.org/resources/openmp-compilers/>


**Compilation**

This is a cmake project. Use the following commands (on Linux) to compile and build the project.

mkdir build

cd build

cmake ..

make

**Example Usage**

./main ../config/params_greenhouse.cfg

The config (.cfg) file passed as an argument to the program can be found in the config/ directory. You can change source and destination directories here, and also tweak parameters like initial mean and standard deviation values for the optimization variables.

Phytomer datasets from simulated, greenhouse and field environments are provided in the data/ directory. When trying new phytomer datasets, make sure to perform the following preprocessing steps on the phytomer point clouds in the dataset,

i. Scale the point cloud so that the phenotype parameter values are in cm units. You can refer to point clouds in the current datasets to get an idea of the scale.

ii. 3D transform the point cloud to be in an approximately similar orientation as the generated phytomers. This is important since Iterative Closest Point (ICP) like algorithms can only locally optimize for the 3D transform. Again you can refer to point clouds in the current datasets to get an idea of this initial transform. 

iii. Downsample the phytomer point clouds to around ~5000 points to limit computation times.

### Who do I talk to? ###

Paloma Sodhi

psodhi@andrew.cmu.edu
