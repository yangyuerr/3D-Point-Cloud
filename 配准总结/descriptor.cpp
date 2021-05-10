/*****************************************************************************
*                                descriptor.cpp
*该文件主要是简单的点云特征提取:
*PFH、FPFH、RoPS、SI、SHOT
*USC、SDASS、SDASS_KP
*
*此文件没有主函数，可在registration.cpp中调用此文件中的函数
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/

//计算LRF
void CAlgorithmSet::rcs_lrf_Z_axis(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, Vertex& z_axis)
{
	int i;
	pcl::PointXYZ query_point = cloud->points[0];

	// 计算协方差矩阵
	Eigen::Matrix3f Cov;
	Eigen::Vector4f centroid;
	pcl::compute3DCentroid(*cloud, centroid);
	pcl::computeCovarianceMatrix(*cloud, centroid, Cov);
	EIGEN_ALIGN16 Eigen::Vector3f::Scalar eigen_min;
	EIGEN_ALIGN16 Eigen::Vector3f normal;
	pcl::eigen33(Cov, eigen_min, normal);
	z_axis.x = normal(0);
	z_axis.y = normal(1);
	z_axis.z = normal(2);
	// 消除z轴方向二义性
	float z_sign = 0;
	for (i = 0; i < cloud->points.size(); i++)
	{
		float vec_x = query_point.x - cloud->points[i].x;
		float vec_y = query_point.y - cloud->points[i].y;
		float vec_z = query_point.z - cloud->points[i].z;
		z_sign += (vec_x * z_axis.x + vec_y * z_axis.y + vec_z * z_axis.z);
	}
	if (z_sign < 0)//正视凸平面
	{
		z_axis.x = -z_axis.x;
		z_axis.y = -z_axis.y;
		z_axis.z = -z_axis.z;
	}
}

//计算法向量
void computeNormals(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_ptr_search, pcl::PointCloud<pcl::Normal>::Ptr output_normals_ptr, double radius)
{
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal>ne;
	ne.setInputCloud(input_cloud_ptr);
	//ne.setNumberOfThreads(4);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree<pcl::PointXYZ>());
	tree_src->setInputCloud(input_cloud_ptr_search);
	ne.setSearchMethod(tree_src);
	ne.setRadiusSearch(radius);
	ne.compute(*output_normals_ptr);

}

//PFH
void computePFH(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::PFHSignature125>::Ptr output_PFH_ptr, double radius)
{
	pcl::PFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::PFHSignature125> pfh_src;
	pfh_src.setInputCloud(input_cloud_ptr);
	pfh_src.setSearchSurface(input_cloud_vol);
	pfh_src.setInputNormals(input_normals_ptr);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
	pfh_src.setSearchMethod(tree_src_fpfh);
	pfh_src.setRadiusSearch(radius);
	pfh_src.compute(*output_PFH_ptr);
	std::cout << "compute *cloud_src pfh" << endl;
}

//FPFH
void computeFPFH(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::FPFHSignature33>::Ptr output_FPFH_ptr, double radius)
{
	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
	fpfh_src.setInputCloud(input_cloud_ptr);
	fpfh_src.setSearchSurface(input_cloud_vol);
	fpfh_src.setInputNormals(input_normals_ptr);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
	fpfh_src.setSearchMethod(tree_fpfh);
	fpfh_src.setRadiusSearch(radius);
	fpfh_src.compute(*output_FPFH_ptr);
}

//RoPS
void computeRoPS(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::VFHSignature308>::Ptr output_RoPS_ptr, double searchRadius, double supportRadius)
{
	output_RoPS_ptr->clear();
	pcl::PointCloud<pcl::PointNormal>::Ptr cloudNormals(new pcl::PointCloud<pcl::PointNormal>);
	pcl::concatenateFields(*input_cloud_vol, *input_normals_ptr, *cloudNormals);
	pcl::search::KdTree<pcl::PointNormal>::Ptr kdtree(new pcl::search::KdTree<pcl::PointNormal>);
	kdtree->setInputCloud(cloudNormals);
	pcl::GreedyProjectionTriangulation<pcl::PointNormal> triangulation;
	pcl::PolygonMesh triangles;
	triangulation.setSearchRadius(0.2);
	triangulation.setMu(2.5);
	triangulation.setMaximumNearestNeighbors(100);
	triangulation.setMaximumSurfaceAngle(M_PI / 4); // 45 degrees.
	triangulation.setNormalConsistency(false);
	triangulation.setMinimumAngle(M_PI / 18); // 10 degrees.
	triangulation.setMaximumAngle(2 * M_PI / 3); // 120 degrees.
	triangulation.setInputCloud(cloudNormals);
	triangulation.setSearchMethod(kdtree);
	triangulation.reconstruct(triangles);

	pcl::ROPSEstimation <pcl::PointXYZ, pcl::Histogram <135> > feature_estimator2;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<pcl::PointXYZ>);
	feature_estimator2.setSearchMethod(tree_tgt_fpfh);
	feature_estimator2.setInputCloud(input_cloud_ptr);
	feature_estimator2.setSearchSurface(input_cloud_vol);
	feature_estimator2.setRadiusSearch(searchRadius);
	feature_estimator2.setNumberOfPartitionBins(5);
	feature_estimator2.setNumberOfRotations(3);
	feature_estimator2.setSupportRadius(supportRadius);
	feature_estimator2.setTriangles(triangles.polygons);

	pcl::PointCloud<pcl::Histogram <135>>::Ptr rops(new pcl::PointCloud <pcl::Histogram <135>>());
	feature_estimator2.compute(*rops);
	pcl::VFHSignature308 midpoint1;
	for (int j = 0; j<rops->size(); ++j)
	{
		for (int i = 0; i<135; i++)
		{
			midpoint1.histogram[i] = rops->points[j].histogram[i];
		}
		for (int i = 135; i < 308; i++)
		{
			midpoint1.histogram[i] = 0;
		}
		output_RoPS_ptr->push_back(midpoint1);
	}
}


//SI
void computeSpinImage(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::VFHSignature308>::Ptr output_SI_ptr, double radius)
{
	output_SI_ptr->clear();
	pcl::SpinImageEstimation<pcl::PointXYZ, pcl::Normal, pcl::Histogram<153>> SI1;
	SI1.setInputCloud(input_cloud_ptr);
	SI1.setSearchSurface(input_cloud_vol);

	SI1.setRadiusSearch(radius);
	SI1.setInputNormals(input_normals_ptr);
	SI1.useNormalsAsRotationAxis();
	SI1.setSupportAngle(0);
	SI1.setRadialStructure(0);
	pcl::PointCloud<pcl::Histogram<153>>::Ptr spinImage(new pcl::PointCloud<pcl::Histogram<153>>());
	SI1.compute(*spinImage);

	pcl::VFHSignature308 midpoint1;
	for (int j = 0; j<spinImage->size(); ++j)
	{
		for (int i = 0; i<153; i++)
		{
			midpoint1.histogram[i] = spinImage->points[j].histogram[i];
		}
		for (int i = 153; i<308; i++)
		{
			midpoint1.histogram[i] =0;
		}
		output_SI_ptr->push_back(midpoint1);
	}
}

//SHOT
void computeSHOT(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::SHOT352>::Ptr output_SHOT_ptr, double radius)
{
	pcl::SHOTEstimation<pcl::PointXYZ, pcl::Normal, pcl::SHOT352> shot_tgt;
	shot_tgt.setInputCloud(input_cloud_ptr);
	shot_tgt.setSearchSurface(input_cloud_vol);
	shot_tgt.setInputNormals(input_normals_ptr);
	//shot_tgt.setNumberOfThreads(4);
	//pcl::PointCloud<pcl::SHOT352>::Ptr fpfhs_tgt(new pcl::PointCloud<pcl::SHOT352>());
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt_shot(new pcl::search::KdTree<pcl::PointXYZ>);
	shot_tgt.setSearchMethod(tree_tgt_shot);
	shot_tgt.setRadiusSearch(radius);
	shot_tgt.compute(*output_SHOT_ptr);

}

//USC
void computeUSC(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr input_cloud_vol, pcl::PointCloud<pcl::Normal>::Ptr input_normals_ptr, pcl::PointCloud<pcl::UniqueShapeContext1960>::Ptr output_USC_ptr, double radius)
{
	pcl::UniqueShapeContext<pcl::PointXYZ, pcl::UniqueShapeContext1960, pcl::ReferenceFrame > usc_tgt;
	usc_tgt.setInputCloud(input_cloud_ptr);
	usc_tgt.setSearchSurface(input_cloud_vol);
	//pcl::PointCloud<pcl::UniqueShapeContext1960>::Ptr fpfhs_tgt(new pcl::PointCloud<pcl::UniqueShapeContext1960>()); // 索引点集对应的usc描述子
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_usc_tgt(new pcl::search::KdTree<pcl::PointXYZ>);
	tree_usc_tgt->setInputCloud(input_cloud_vol);
	usc_tgt.setSearchMethod(tree_usc_tgt);
	usc_tgt.setRadiusSearch(radius);
	usc_tgt.setMinimalRadius(radius / 10.0);
	// Radius used to compute the local point density for the neighbors
	// (the density is the number of points within that radius).
	usc_tgt.setPointDensityRadius(radius / 5.0);
	usc_tgt.compute(*output_USC_ptr);
}

void SDASS(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr& keypoint, float& feat_radius, float& lma_sup_radius, int depth_bins, int radius_bins, int lraI_bins, pcl::PointCloud<pcl::ESFSignature640>::Ptr SDASS)
{
	SDASS->clear();
	int i, j, m;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	vector<int>pointIdx, temp1;
	vector<float>pointDst, temp2;
	kdtree.setInputCloud(cloud);
	vector<vector<float>> histogram;
	// (1)For all points compute the unambiguous LRA, and use the LRA to compute the Normal included angle
	vector<Vertex> LRA_of_all_points;
	for (i = 0; i < cloud->points.size(); i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr pt_lra(new pcl::PointCloud<pcl::PointXYZ>);
		Vertex lra;
		kdtree.radiusSearch(cloud->points[i], lma_sup_radius, pointIdx, pointDst);
		for (j = 0; j < pointIdx.size(); j++)
		{
			pt_lra->points.push_back(cloud->points[pointIdx[j]]);//确定局部曲面
		}
		rcs_lrf_Z_axis(pt_lra, lra);
		LRA_of_all_points.push_back(lra);
	}

	// (2)For all keypoints, compute the SDASS descriptor

	for (i = 0; i < keypoint->size(); i++)
	{
		vector<float> histogram_kp;
		vector<int>neiIdx;
		vector<float>neiDst;
		Vertex lra;
		pcl::PointCloud<pcl::PointXYZ>::Ptr pt_lra(new pcl::PointCloud<pcl::PointXYZ>);
		kdtree.radiusSearch(keypoint->points[i], lma_sup_radius, pointIdx, pointDst);
		for (j = 0; j < pointIdx.size(); j++)
		{
			pt_lra->points.push_back(cloud->points[pointIdx[j]]);//确定局部曲面
		}
		rcs_lrf_Z_axis(pt_lra, lra);

		//search the neighbor points
		kdtree.radiusSearch(keypoint->points[i], feat_radius, neiIdx, neiDst);
		SDASS_KP(cloud, feat_radius, keypoint->points[i], lra, neiIdx, LRA_of_all_points, depth_bins, radius_bins, lraI_bins, histogram_kp);
		histogram.push_back(histogram_kp);
	}

	//pcl::PointCloud<pcl::ESFSignature640>::Ptr SDASS(new pcl::PointCloud<pcl::ESFSignature640>);
	pcl::ESFSignature640 midpoint;

	for (int j = 0; j < histogram.size(); ++j)
	{
		for (int i = 0; i < histogram[j].size(); i++)
		{
			midpoint.histogram[i] = histogram[j][i];
		}
		for (int i = histogram[j].size(); i < 512; i++)
		{
			midpoint.histogram[i] = 0;

		}
		SDASS->push_back(midpoint);
	}
}

void SDASS_KP(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, float& radius, pcl::PointXYZ& kp, Vertex kp_lra, vector<int>& neigh_Idx,
	vector<Vertex>& LRAs, int& num_depth_bins, int& num_radius_bins, int& num_normalangle_bins, vector<float>& histogram)
{
	int i;
	float nx = kp_lra.x, ny = kp_lra.y, nz = kp_lra.z;   // lra of kpt
	float x0 = kp.x - nx * radius, y0 = kp.y - ny * radius, z0 = kp.z - nz * radius;

	vector<float> SDASS_histogram(num_depth_bins * num_radius_bins * num_normalangle_bins, 0);

	float plane_D = -(nx * x0 + ny * y0 + nz * z0);
	int depth_bin_num = num_depth_bins;
	float depth_stride = 2 * radius / depth_bin_num;
	int depth_bin_id;

	//////////////////////////////////////////////////////////////////////////
	int density_bin_num = num_radius_bins;
	float density_stride = radius / density_bin_num;
	int density_bin_id;

	//////////////////////////////////////////////////////////////////////////
	int angle_bin_num = num_normalangle_bins;
	float angle_stride = 180.0f / angle_bin_num;
	int angle_bin_id;

	//////////////////////////////////////////////////////////////////////////
	float a, b, c;
	float temp_depth, temp_radius, temp_angle;
	for (i = 0; i < neigh_Idx.size(); i++)
	{

		pcl::PointXYZ neighPt = cloud->points[neigh_Idx[i]];
		temp_depth = nx * neighPt.x + ny * neighPt.y + nz * neighPt.z + plane_D;
		c = (neighPt.x - kp.x) * (neighPt.x - kp.x) +
			(neighPt.y - kp.y) * (neighPt.y - kp.y) +
			(neighPt.z - kp.y) * (neighPt.z - kp.y);
		b = (neighPt.x - kp.x) * nx + (neighPt.y - kp.y) * ny + (neighPt.z - kp.z) * nz;
		a = sqrt(abs(c - b * b));

		temp_radius = a;
		temp_angle = LRAs[neigh_Idx[i]].x * nx + LRAs[neigh_Idx[i]].y * ny + LRAs[neigh_Idx[i]].z * nz;
		if (temp_angle > 1) temp_angle = 1;
		if (temp_angle < -1) temp_angle = -1;
		temp_angle = acos(temp_angle) / M_PI * 180;

		//compute histograms
		float depth_bin_temp = temp_depth / depth_stride;
		int depth_bin = depth_bin_temp + 1;
		if (depth_bin < 1) depth_bin = 1;
		if (depth_bin > depth_bin_num) depth_bin = depth_bin_num;

		float radius_bin_temp = temp_radius / density_stride;
		int radius_bin = radius_bin_temp + 1;
		if (radius_bin < 1) radius_bin = 1;
		if (radius_bin > density_bin_num) radius_bin = density_bin_num;

		float angle_bin_temp = temp_angle / angle_stride;
		int angle_bin = angle_bin_temp + 1;
		if (angle_bin < 1) angle_bin = 1;
		if (angle_bin > angle_bin_num) angle_bin = angle_bin_num;

		SDASS_histogram[(depth_bin - 1) * num_radius_bins * num_normalangle_bins + (radius_bin - 1) * num_normalangle_bins + angle_bin - 1] += 1 / float(neigh_Idx.size());
	}
	histogram = SDASS_histogram;
}
