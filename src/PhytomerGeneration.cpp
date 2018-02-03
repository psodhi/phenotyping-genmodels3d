#include "PhytomerGeneration.h"

PhytomerGeneration::PhytomerGeneration()
{
	lsys_string_ = "F[+G]FX"; 	// string hardcoded for a phytomer (stage_number = 0)
	scale_ = 1;
	leaf_area_ = 0;

	sample_points_ = 6000;
	percentage_stem_points_ = 0.4;
}

void PhytomerGeneration::updateParameters(float stem_radius_arg, float leaf_angle_arg, float leaf_length_arg)
{
	// 3 Variables in CE Optimization set externally
	stem_radius_ = stem_radius_arg;
	leaf_angle_ = leaf_angle_arg;
	leaf_length_ = leaf_length_arg;

	leaf_width_ = (leaf_length_/5)* 0.8;
	leaf_rotate_cont_angle_ = 3;
	stem_internode_length_ = leaf_length_/4 + 1.5;

	leaf_thickness_ = 1;
}

void PhytomerGeneration::updateParameters(float stem_radius_arg, float leaf_angle_arg, float leaf_length_arg, float stem_internode_length_arg, float leaf_rotate_cont_angle_arg, float leaf_width_arg)
{
	// All 6 Variables in CE Optimization set externally
	stem_radius_ = stem_radius_arg;
	leaf_angle_ = leaf_angle_arg;
	leaf_length_ = leaf_length_arg;
	stem_internode_length_ = stem_internode_length_arg;
	leaf_rotate_cont_angle_ = leaf_rotate_cont_angle_arg;
	leaf_width_ = leaf_width_arg;

	// TODO : add these variables to Cross Entropy Optimization
	leaf_thickness_ = 1;
}

void PhytomerGeneration::setNumberOfSamplePoints(int sample_points_arg)
{
	sample_points_ = sample_points_arg;
}

pcl::PointCloud<PointT> PhytomerGeneration::getPhytomerPointCloud()
{
	return cloud_gen_phytomer_;
}

// Return Noisy Point Cloud (with std_dev_perc*100% gaussian noise added to each 3D point)
pcl::PointCloud<PointT> PhytomerGeneration::getPhytomerPointCloudNoisy(float std_dev_perc)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);

	for (int i = 0; i < cloud_gen_phytomer_.points.size(); i++)
	{
		float v = pow( pow(cloud_gen_phytomer_.points[i].x, 2)+pow(cloud_gen_phytomer_.points[i].y, 2)+pow(cloud_gen_phytomer_.points[i].z, 2) , 0.5); 

		std::normal_distribution<double> distribution_x (cloud_gen_phytomer_.points[i].x, std_dev_perc*cloud_gen_phytomer_.points[i].x);
		// std::normal_distribution<double> distribution_x (cloud_gen_phytomer_.points[i].x, std_dev_perc*v);
		cloud_gen_phytomer_.points[i].x = distribution_x(generator);

		std::normal_distribution<double> distribution_y (cloud_gen_phytomer_.points[i].y, std_dev_perc*cloud_gen_phytomer_.points[i].y);
		// std::normal_distribution<double> distribution_y (cloud_gen_phytomer_.points[i].y, std_dev_perc*v);
		cloud_gen_phytomer_.points[i].y = distribution_y(generator);

		std::normal_distribution<double> distribution_z (cloud_gen_phytomer_.points[i].z, std_dev_perc*cloud_gen_phytomer_.points[i].z);
		// std::normal_distribution<double> distribution_z (cloud_gen_phytomer_.points[i].z, std_dev_perc*v);
		cloud_gen_phytomer_.points[i].z = distribution_z(generator);
	}

	return cloud_gen_phytomer_;
}

PhytomerGeneration::~PhytomerGeneration()
{
}
Eigen::Matrix3f PhytomerGeneration::make_rotation_matrix(float leaf_angle_rad, char symbol)
{
	Eigen::Matrix3f rotation_matrix;

	switch(symbol)
	{
		case '+': // anticlockwise about z
			rotation_matrix << cos(leaf_angle_rad), sin(leaf_angle_rad), 0,
						-sin(leaf_angle_rad), cos(leaf_angle_rad), 0,
						0, 0, 1;
		break;

		case '-': // clockwise about z
			rotation_matrix << cos(-leaf_angle_rad), sin(-leaf_angle_rad), 0,
						-sin(-leaf_angle_rad), cos(-leaf_angle_rad), 0,
						0, 0, 1;
		break;

		default:
			cout << "ERROR: Invalid symbol passed to make_rotation_matrix()" << endl;
	}

	return rotation_matrix;
}

float PhytomerGeneration::deg2rad(float angle_degrees)
{
	float angle_radians = angle_degrees*(M_PI/180);
	return angle_radians;
}

float PhytomerGeneration::fluctuate(float x, float ratio)
{
	// Given x and ratio, generate a number between [x(1-ratio), x(1+ratio)]
	float x_fluctuate = (rand() * 2 * x * ratio - (x * ratio)) + x;
	return x_fluctuate;
}

tuple<Eigen::Matrix<float, 4, 3> , Eigen::Vector3f> PhytomerGeneration::get_leaf_cross_section(Eigen::Vector3f pt, int idx, Eigen::Vector3f piece_normal, Eigen::Vector3f cp_backup, Eigen::Vector3f v, int pieces_per_leaf)
{
	// Compute Initial Variables
	Eigen::Vector3f z_vec(0, 0.15 + 0.7, 0.3 - 0.3);
	z_vec = z_vec / z_vec.norm();

	float concavity = 2; // [1,3]
	// Control these two vars for fixing the polygonal patches near end of leaf
    float leaf_first_max_width_point = floor(fluctuate(pieces_per_leaf / 5, 0)) + 5; // [8, 13]
    float leaf_last_max_width_point = pieces_per_leaf - floor(fluctuate(pieces_per_leaf / 5, 0)) - 5; // [10, 20]

    Eigen::Matrix<float, 4, 3> spts;
    Eigen::Vector3f cp;

    spts(0,0) = pt(0); spts(0,1) = pt(1); spts(0,2) = pt(2);
    spts(1,0) = pt(0)-leaf_thickness_*piece_normal(0); spts(1,1) = pt(1)-leaf_thickness_*piece_normal(1); spts(1,2) = pt(2)-leaf_thickness_*piece_normal(2);

    cp = z_vec.cross(v);
    cp = cp / cp.norm();

    if (v(1) < -0.8)
        cp = cp_backup;

    cp(0) = 0; cp(1) = 0; cp(2) = 1; 

    float part_width;
    if (idx < leaf_first_max_width_point)
    {
    	// i. concave_increasing_half function
    	// float max_value = leaf_width_ / 2;
    	// float ratio =  idx/(pieces_per_leaf/3);
    	// float k = 2 + concavity;
    	// part_width = 1/(1 - exp(-k) + exp(-k * ratio)) * max_value;

    	// ii. fixed_width_half function
    	part_width = leaf_width_/2;
    }
    else if (idx >= leaf_first_max_width_point)
    {
    	// concave_decreasing function
    	float max_value = leaf_width_ / 2;
    	float ratio = ((idx - leaf_last_max_width_point) / (pieces_per_leaf - leaf_last_max_width_point));
    	float k = concavity;
    	float local_max = (k + log(1 + 1/exp(k)))/k;
    	part_width = ((k + log(1-ratio + (1/exp(k))))/k / local_max) * max_value;
    	// part_width = leaf_width_/2;
    }
    else
    	part_width = leaf_width_/2;

    spts(2, 0) = pt(0) + part_width * cp(0); spts(2, 1) = pt(1) + part_width * cp(1); spts(2, 2) = pt(2) + part_width * cp(2);
    spts(3, 0) = pt(0) + part_width * (-cp(0)); spts(3, 1) = pt(1) + part_width * (-cp(1)); spts(3, 2) = pt(2) + part_width * (-cp(2));

    tuple<Eigen::Matrix<float, 4, 3> , Eigen::Vector3f> tuple_temp = make_tuple(spts, cp);
    return tuple_temp;
}

void PhytomerGeneration::make_leaf_patches(vector<endverts> vertsG_new, int pieces_per_leaf, vtkSmartPointer<vtkPolyData> &polydata_leaf)
{
	int number_leaves = vertsG_new.size() / pieces_per_leaf;
	Eigen::Vector3f vz(0.0, 1.0, 0.0);

	for (int i = 0; i < number_leaves; i++)
	{
		Eigen::Vector3f v;
		Eigen::Vector3f piece_normal;
		Eigen::Vector3f prev_normal(0.0, 1.0, 0.0);
		Eigen::Vector3f cp_backup(0.0, 0.0, 1.0); // originally, cp_backup was [0 0 1]

		bool leaf_broken;
		int cp_idx = 14;

		vtkSmartPointer<vtkPoints> points_leaf = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> polygons_leaf = vtkSmartPointer<vtkCellArray>::New();

		for (int j = 0; j < vertsG_new.size(); j++)
		{
			v = vertsG_new[j].endpos - vertsG_new[j].startpos;
			v = v / v.norm();

			piece_normal = -(v + vz * (-v.dot(v)/v.dot(vz)));
			piece_normal = piece_normal / piece_normal.norm();
			float normal_angle = acos(vz.dot(piece_normal));

			// Make sure it is pointing downwards
			if (normal_angle > M_PI/2)
				piece_normal = piece_normal; // DEBUG, originally : piece_normal = piece_normal

			if ( (normal_angle > M_PI/2) || acos(prev_normal.dot(piece_normal)) > (M_PI/2) )
				leaf_broken = true;

			// Get set of 4 points for both start_section and end_section
			tuple<Eigen::Matrix<float, 4, 3> , Eigen::Vector3f> tuple_temp;
			tuple_temp = get_leaf_cross_section(vertsG_new[j].startpos, j+1, piece_normal, cp_backup, v, pieces_per_leaf);
			Eigen::Matrix<float, 4, 3> start_section = get<0>(tuple_temp); Eigen::Vector3f cp = get<1>(tuple_temp);
			tuple_temp = get_leaf_cross_section(vertsG_new[j].endpos, j+1, piece_normal, cp_backup, v, pieces_per_leaf);
			Eigen::Matrix<float, 4, 3> end_section = get<0>(tuple_temp);

			if (j+1 == cp_idx)
				cp_backup = cp;

			// quad surf 1
			vtkSmartPointer<vtkPolygon> polygon1 = vtkSmartPointer<vtkPolygon>::New();

			points_leaf->InsertNextPoint(start_section(0,0), start_section(0,1), start_section(0,2));
			points_leaf->InsertNextPoint(end_section(0,0), end_section(0,1), end_section(0,2));
			points_leaf->InsertNextPoint(end_section(3,0), end_section(3,1), end_section(3,2));
			points_leaf->InsertNextPoint(start_section(3,0), start_section(3,1), start_section(3,2));
			polygon1->GetPointIds()->SetNumberOfIds(4); //make a quad
			polygon1->GetPointIds()->SetId(0, 8*j+0);
			polygon1->GetPointIds()->SetId(1, 8*j+1);
			polygon1->GetPointIds()->SetId(2, 8*j+2);
			polygon1->GetPointIds()->SetId(3, 8*j+3);

			polygons_leaf->InsertNextCell(polygon1);

			// quad surf 2
			vtkSmartPointer<vtkPolygon> polygon2 = vtkSmartPointer<vtkPolygon>::New();

			points_leaf->InsertNextPoint(start_section(0,0), start_section(0,1), start_section(0,2));
			points_leaf->InsertNextPoint(end_section(0,0), end_section(0,1), end_section(0,2));
			points_leaf->InsertNextPoint(end_section(2,0), end_section(2,1), end_section(2,2));
			points_leaf->InsertNextPoint(start_section(2,0), start_section(2,1), start_section(2,2));
			polygon2->GetPointIds()->SetNumberOfIds(4); //make a quad
			polygon2->GetPointIds()->SetId(0, 8*j+4);
			polygon2->GetPointIds()->SetId(1, 8*j+5);
			polygon2->GetPointIds()->SetId(2, 8*j+6);
			polygon2->GetPointIds()->SetId(3, 8*j+7);

			polygons_leaf->InsertNextCell(polygon2);

		}
		polydata_leaf->SetPoints(points_leaf);
		polydata_leaf->SetPolys(polygons_leaf);
	}
	// visualizeVTKPolyData(polydata_leaf);
}

void PhytomerGeneration::plot_leaves(vector<endverts> vertsG, vtkSmartPointer<vtkPolyData> &polydata_leaf, pcl::PointCloud<PointT>::Ptr cloud_out_leaf)
{
	int num_vertices = vertsG.size();
	int pieces_per_leaf = 30;
	float dot_epsilon = 0.9; // controls when to change rotation rate of leaf
	Eigen::Vector3f h_vec(0.0, 1.0, 0.0);

	float rotate_cont_angle, rotate_small_angle;
	rotate_cont_angle = this->fluctuate(deg2rad(leaf_rotate_cont_angle_), 0.0) + deg2rad(1.5);
	rotate_small_angle = this->fluctuate(deg2rad(0.5*leaf_rotate_cont_angle_), 0.0);

	endverts vert;
	vector<endverts> vertsG_new;
	Eigen::Vector3f new_pos;
	Eigen::Vector3f start_pos = vertsG[0].startpos; 

	for (int i = 0; i < vertsG.size() ; i++)
	{
		float step_length = leaf_length_ / pieces_per_leaf;
		Eigen::Vector3f v = vertsG[0].endpos - vertsG[0].startpos;
		v = v / v.norm();

		for (int j = 0; j < pieces_per_leaf; j++)
		{
			// accumulate new vector of vertices composed of smaller line segments
			new_pos = start_pos + v * step_length;
			vert.startpos = start_pos; vert.endpos = new_pos;
			vertsG_new.push_back(vert);
			start_pos = new_pos;

			// compute rotation vector for bending vector v
			// when the leaf is close to vertical line, reduce the bending rate
			Eigen::Quaternion<float> rotvec_quat = Eigen::Quaternion<float>::FromTwoVectors(h_vec, v);

			Eigen::Vector3f rotate_axis(rotvec_quat.x(), rotvec_quat.y(), rotvec_quat.z());
			rotate_axis = rotate_axis / rotate_axis.norm();

			Eigen::AngleAxis<float> rotvec_angleaxis;
			if ( v.dot(h_vec) >= dot_epsilon )
				rotvec_angleaxis = Eigen::AngleAxis<float>(rotate_small_angle, rotate_axis);
			else
				rotvec_angleaxis = Eigen::AngleAxis<float>(rotate_cont_angle, rotate_axis);

	        // update vector v to bend the leaf
			v = rotvec_angleaxis.toRotationMatrix() * v;
			v = v / v.norm();
		}

	}

	this->make_leaf_patches(vertsG_new, pieces_per_leaf, polydata_leaf);

	// Convert surface mesh to triangle mesh representation
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New ();
	triangleFilter->SetInputData(polydata_leaf);
	triangleFilter->Update();
	polydata_leaf = triangleFilter->GetOutput();

	pcl::PointCloud<PointT>::Ptr cloud_mesh_sampled (new pcl::PointCloud<PointT>);
	PointsFromSurfaceSampler sampler1;
	sampler1.uniform_sampling(polydata_leaf, (1-percentage_stem_points_)*sample_points_, 0, *cloud_mesh_sampled);

	*cloud_out_leaf = *cloud_out_leaf + *cloud_mesh_sampled;
}

void PhytomerGeneration::plot_stem_nodes(vector<endverts> vertsF, vtkSmartPointer<vtkPolyData> &polydata_stem, pcl::PointCloud<PointT>::Ptr cloud_out_stem)
{
	for (int i = 0; i < vertsF.size(); i++)
	{
		Eigen::Vector3f centerpos = 0.5*(vertsF[i].startpos + vertsF[i].endpos);
		vtkSmartPointer<vtkCylinderSource> cylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
		cylinderSource->SetCenter(centerpos(0), centerpos(1), centerpos(2));
		cylinderSource->SetRadius(stem_radius_);
		cylinderSource->SetHeight(stem_internode_length_);
		cylinderSource->SetResolution(100);
		cylinderSource->CappingOff();

		// Convert surface mesh to triangle mesh representation
		vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New ();
		triangleFilter->SetInputConnection(cylinderSource->GetOutputPort());
		triangleFilter->Update();
		polydata_stem = triangleFilter->GetOutput();

		pcl::PointCloud<PointT>::Ptr cloud_mesh_sampled (new pcl::PointCloud<PointT>);
		PointsFromSurfaceSampler sampler1;
		sampler1.uniform_sampling(polydata_stem, percentage_stem_points_*sample_points_, 0, *cloud_mesh_sampled);

		*cloud_out_stem = *cloud_out_stem + *cloud_mesh_sampled;
	}
}

void PhytomerGeneration::plot_from_string()
{
	float leaf_angle_rad = deg2rad(leaf_angle_);

	Eigen::Vector3f position(0.0, 0.0, 0.0);
	Eigen::Vector3f orientation(0.0, 1.0, 0.0);
	Eigen::Vector3f new_position;
	Eigen::Matrix3f rotation_matrix;

	// stem, leaf composed of vector of line segments
	// vertsF, vertsG hold start and end vertex positions for each line segment making up stem, leaf respectively  
	vector<endverts> vertsF, vertsG;
	endverts vert;
	Eigen::Vector3f pos;
	Eigen::Vector3f ori_proj(1.0, 0.0, 0.0);

	pose pose_curr;
	stack<pose> pose_stack;

	for (int i = 0; i < lsys_string_.size(); i++)
	{
		switch(lsys_string_[i])
		{
			case 'F' : // collect stem vertices
				new_position = position + stem_internode_length_*orientation;
				vert.startpos = position; vert.endpos = new_position;
				vertsF.push_back(vert);

				position = new_position;
			break;

			case 'G' : // collect leaf vertices
				pos = position + stem_radius_ * ori_proj;
				new_position = pos + leaf_length_ * (orientation/orientation.norm());

				vert.startpos = pos; vert.endpos = new_position;
				vertsG.push_back(vert);

				position = new_position;
			break;

			case '+' :
				rotation_matrix = make_rotation_matrix(leaf_angle_rad, '+');
				orientation = rotation_matrix * orientation;
			break;

			case '[' :
				// save the current pose and push on the stack
				pose_curr.position = position;
				pose_curr.orientation = orientation;
				pose_stack.push(pose_curr);
			break;

			case ']' :
				// pop on the stack and return to the saved position
				pose_curr = pose_stack.top();
				position = pose_curr.position; 
				orientation = pose_curr.orientation;
				pose_stack.pop();
			break;

			case 'X' :
				// placeholder character : do nothing
			break;

			default :
				cout << "ERROR : Invalid characters in L-system model string" << endl;
			break;			
		}
	}

	pcl::PointCloud<PointT>::Ptr cloud_stem	(new pcl::PointCloud<PointT>);
	pcl::PointCloud<PointT>::Ptr cloud_leaf	(new pcl::PointCloud<PointT>);

	vtkSmartPointer<vtkPolyData> polydata_stem = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> polydata_leaf = vtkSmartPointer<vtkPolyData>::New();

	plot_stem_nodes(vertsF, polydata_stem, cloud_stem);
	plot_leaves(vertsG, polydata_leaf, cloud_leaf);

	// Colorize cloud_stem and cloud_leaf differently before accumulating into cloud_gen_phytomer_
	int rgb_vals_stem[3] = {100, 0, 0}; 
	colorizeCloud(*cloud_stem, rgb_vals_stem);
 	int rgb_vals_leaf[3] = {0, 100, 0};
	colorizeCloud(*cloud_leaf, rgb_vals_leaf);

	cloud_gen_phytomer_.clear();
	cloud_gen_phytomer_ = *cloud_stem + *cloud_leaf;
	cloud_gen_phytomer_.width = cloud_gen_phytomer_.points.size(); 
	cloud_gen_phytomer_.height = 1;

	Eigen::Vector4f pt_centroid;
	compute3DCentroid (cloud_gen_phytomer_, pt_centroid);
	demeanPointCloud (cloud_gen_phytomer_, pt_centroid, cloud_gen_phytomer_);

}

void PhytomerGeneration::generate_sorghum_phytomer()
{
	plot_from_string();
}


void PhytomerGeneration::colorizeCloud(pcl::PointCloud<PointT> &cloud_in, int rgb_vec[])
{
	for (int i = 0; i < cloud_in.points.size(); i++)
	{
		cloud_in.points[i].r = rgb_vec[0];
		cloud_in.points[i].g = rgb_vec[1];
		cloud_in.points[i].b = rgb_vec[2];
	}
}

int PhytomerGeneration::visualizeVTKPolyData(vtkSmartPointer<vtkPolyData> polydata_vis)
{
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputDataObject(polydata_vis);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	//Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	
	// Add the actor to the scene
	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1); // Background color dark green
	// renderer->SetBackground(.1, .3,.2); // Background color dark green
	
	// Render and interact
	renderWindow->SetWindowName("polydata_vis");
	renderWindow->Render();
	renderWindowInteractor->Start();
	
	return EXIT_SUCCESS;
}

int PhytomerGeneration::visualizeVTKPolyData(vtkSmartPointer<vtkPolyData> polydata_vis_1, vtkSmartPointer<vtkPolyData> polydata_vis_2)
{
	// initialize actor1 with polydata_vis_1
	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputDataObject(polydata_vis_1);
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper1);
	actor1->GetProperty()->SetColor(0.6, 0, 0);

	// initialize actor2 with polydata_vis_2
	vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputDataObject(polydata_vis_2);
	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetColor(0, 0.4, 0);

	//Create a renderer, render window, and interactor
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	
	// Add the actor to the scene
	renderer->AddActor(actor1);
	renderer->AddActor(actor2);
	renderer->SetBackground(1, 1, 1); // Background color dark green
	// renderer->SetBackground(.1, .3,.2); // Background color dark green
	
	// Render and interact
	renderWindow->SetWindowName("polydata_vis");
	renderWindow->Render();
	renderWindowInteractor->Start();
	
	return EXIT_SUCCESS;
}

int PhytomerGeneration::main()
{

	vtkSmartPointer<vtkCylinderSource> cylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
	cylinderSource->SetCenter(0.0, 0.0, 0.0);
	cylinderSource->SetRadius(stem_radius_);
	cylinderSource->SetHeight(stem_internode_length_);
	cylinderSource->SetResolution(100);

	vtkSmartPointer<vtkPolyData> cylinderPolyData = vtkSmartPointer<vtkPolyData>::New();
	cylinderPolyData = cylinderSource->GetOutput();

	// Convert surface mesh to triangle mesh representation
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New ();
	triangleFilter->SetInputData(cylinderPolyData);
	triangleFilter->Update();
	cylinderPolyData = triangleFilter->GetOutput();

	PointsFromSurfaceSampler sampler1;

	pcl::PointCloud<PointT>::Ptr cloud_mesh_sampled (new pcl::PointCloud<PointT>);
	sampler1.uniform_sampling(cylinderPolyData, 5000, 0, *cloud_mesh_sampled);

	return 0;
}