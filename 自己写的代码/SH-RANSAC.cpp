/*****************************************************************************
*                                SH-RANSAC.cpp
*
*SH-RANSAC是作者在解决实际问题时提出的算法----单点Hough-RANSAC 简称 SH-RANSAC
*首先利用每对关键点的LRF即可得到变换矩阵
*接着对每个关键点进行Hough投票，选出点数最多的格子，
*然后遍历格子中所有关键点对的变换矩阵，变换源点云，并加上ICP算法
*最后根据重合率与RMSE指标选择最优配准结果
*
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/
double applyHeuristicAlign(pcl::PointCloud<PointT>::Ptr cloud_vol, pcl::PointCloud<PointT>::Ptr Keypoint_cloud, pcl::PointCloud<PointT>::Ptr cloud_tgt_vol, pcl::PointCloud<PointT>::Ptr Keypoint_tgt_cloud,
	pcl::PointCloud<pcl::PFHSignature125>::Ptr input_src_feature, pcl::PointCloud<pcl::PFHSignature125>::Ptr input_tgt_feature, pcl::PointCloud<pcl::Normal>::Ptr normals, pcl::PointCloud<pcl::Normal>::Ptr normals_tgt,
	Eigen::Matrix4f& rough_output, Eigen::Matrix4f& fine_output, double &overlap)
{
	pcl::CorrespondencesPtr model_scene_corrs(new pcl::Correspondences());
	pcl::KdTreeFLANN<pcl::PFHSignature125> match_search;
	match_search.setInputCloud(input_tgt_feature);

	//  For each scene keypoint descriptor, find nearest neighbor into the model keypoints descriptor cloud and add it to the correspondences vector.
	for (size_t i = 0; i < input_src_feature->size(); ++i)
	{

		std::vector<int> neigh_indices(2);
		std::vector<float> neigh_sqr_dists(2);
		std::vector<int> neigh_indices_tgt(2);
		std::vector<float> neigh_sqr_dists_tgt(2);
		if (!std::isfinite(input_src_feature->at(i).histogram[0])) //skipping NaNs
		{
			continue;
		}
		int found_neighs = match_search.nearestKSearch(input_src_feature->at(i), 2, neigh_indices, neigh_sqr_dists);
		//pcl::Correspondence corr(static_cast<int> (i), neigh_indices[0], neigh_sqr_dists[0]);
		if (found_neighs == 2 && neigh_sqr_dists[0] / neigh_sqr_dists[1] < 0.8) //  add match only if the squared descriptor distance is less than 0.25 (SHOT descriptor distances are between 0 and 1 by design)
		{
			pcl::Correspondence corr(static_cast<int> (i), neigh_indices[0], neigh_sqr_dists[0]);
			model_scene_corrs->push_back(corr);

		}
	}
	pcl::PointCloud<RFType>::Ptr model_rf(new pcl::PointCloud<RFType>());
	pcl::PointCloud<RFType>::Ptr scene_rf(new pcl::PointCloud<RFType>());

	float rf_rad_(0.35f);
	//特征估计的方法（点云，法线，参考帧）
	pcl::BOARDLocalReferenceFrameEstimation<PointT, NormalType, RFType> rf_est;
	rf_est.setFindHoles(true);
	rf_est.setRadiusSearch(rf_rad_);   //设置搜索半径

	rf_est.setInputCloud(Keypoint_tgt_cloud);  //模型关键点
	rf_est.setInputNormals(normals_tgt); //模型法线
	rf_est.setSearchSurface(cloud_tgt_vol);    //模型
	rf_est.compute(*model_rf);      //模型的参考帧   

	rf_est.setInputCloud(Keypoint_cloud);  //同理
	rf_est.setInputNormals(normals);
	rf_est.setSearchSurface(cloud_vol);
	rf_est.compute(*scene_rf);
	double cg_size_(1.f);
	float cg_thresh_(-0.001f);

	std::vector<Eigen::Matrix4f, Eigen::aligned_allocator<Eigen::Matrix4f> > rototranslations;
	std::vector<pcl::Correspondences> clustered_corrs;

	pcl::Hough3DGrouping<PointT, PointT, RFType, RFType> clusterer;
	clusterer.setHoughBinSize(cg_size_);//霍夫空间设置每个bin的大小
	clusterer.setHoughThreshold(cg_thresh_);//阀值
	clusterer.setUseInterpolation(true);
	clusterer.setUseDistanceWeight(false);

	clusterer.setInputCloud(Keypoint_cloud);
	clusterer.setInputRf(scene_rf);   //设置输入的参考帧
	clusterer.setSceneCloud(Keypoint_tgt_cloud);
	clusterer.setSceneRf(model_rf);
	clusterer.setModelSceneCorrespondences(model_scene_corrs);//model_scene_corrs存储配准的点

	//clusterer.cluster (clustered_corrs);
	//辨认出聚类的对象
	clusterer.recognize(rototranslations, clustered_corrs);


	model_scene_corrs->clear();
	//model_scene_corrs->push_back(clustered_corrs[0][4]);

	for (int i = 0; i < clustered_corrs.size(); i++)
	{
		for (int j = 0; j < clustered_corrs[i].size(); j++)
			model_scene_corrs->push_back(clustered_corrs[i][j]);
	}

	clock_t rf_end = clock();

	//std::cout << "rf time:  " << (rf_end - rf_start)/1000. << endl;

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_rac_trans(new  pcl::PointCloud<pcl::PointXYZ>);
	Eigen::Matrix4f transfor_cor_best = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f transfor_cor_icp = Eigen::Matrix4f::Identity();
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_icp_trans(new  pcl::PointCloud<pcl::PointXYZ>);
	double score_best = 1;
	double calTrans_time = 0;
	double icp_time = 0;
	for (int i = 0; i < model_scene_corrs->size(); i++)
	{
		clock_t calTrans_start = clock();
		pcl::ReferenceFrame scene_point_rf = scene_rf->at(model_scene_corrs->at(i).index_query);
		pcl::ReferenceFrame model_point_rf = model_rf->at(model_scene_corrs->at(i).index_match);

		//去除NAN
		if (isnan(scene_point_rf.x_axis[0]) || isnan(model_point_rf.x_axis[0]))continue;
		//RT矩阵的计算
		Eigen::Matrix3f Rotation_sce;
		Rotation_sce << scene_point_rf.x_axis[0], scene_point_rf.y_axis[0], scene_point_rf.z_axis[0],
			scene_point_rf.x_axis[1], scene_point_rf.y_axis[1], scene_point_rf.z_axis[1],
			scene_point_rf.x_axis[2], scene_point_rf.y_axis[2], scene_point_rf.z_axis[2];//scene_rf->at(model_scene_corrs->at(0).index_query);

		Eigen::Matrix3f Rotation_tgt;
		Rotation_tgt << model_point_rf.x_axis[0], model_point_rf.y_axis[0], model_point_rf.z_axis[0],
			model_point_rf.x_axis[1], model_point_rf.y_axis[1], model_point_rf.z_axis[1],
			model_point_rf.x_axis[2], model_point_rf.y_axis[2], model_point_rf.z_axis[2];
		Eigen::Vector3f tr;
		Eigen::Vector3f tg;
		tg << Keypoint_tgt_cloud->points[model_scene_corrs->at(i).index_match].x, Keypoint_tgt_cloud->points[model_scene_corrs->at(i).index_match].y, Keypoint_tgt_cloud->points[model_scene_corrs->at(i).index_match].z;
		Eigen::Vector3f sr;
		sr << Keypoint_cloud->points[model_scene_corrs->at(i).index_query].x, Keypoint_cloud->points[model_scene_corrs->at(i).index_query].y, Keypoint_cloud->points[model_scene_corrs->at(i).index_query].z;
		Eigen::Matrix3f Rotation_final;
		Rotation_final = Rotation_tgt * Rotation_sce.transpose();
		tr = tg - Rotation_final * sr;
		Eigen::Matrix4f transfor_cor;
		transfor_cor << Rotation_final(0, 0), Rotation_final(0, 1), Rotation_final(0, 2), tr(0, 0),
			Rotation_final(1, 0), Rotation_final(1, 1), Rotation_final(1, 2), tr(1, 0),
			Rotation_final(2, 0), Rotation_final(2, 1), Rotation_final(2, 2), tr(2, 0),
			0, 0, 0, 1;

		//icp配准
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_final(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
		icp.setInputSource(cloud_vol);
		icp.setInputTarget(cloud_tgt_vol);
		// 最大迭代次数
		icp.setMaximumIterations(50);
		// 两次变化矩阵之间的差值
		icp.setTransformationEpsilon(1e-10);
		// 均方误差
		icp.setEuclideanFitnessEpsilon(0.02);
		icp.align(*cloud_final, transfor_cor);
		double score = icp.getFitnessScore();
		clock_t icp_end = clock();
		if (score_best > score)
		{
			score_best = score;
			transfor_cor_icp = icp.getFinalTransformation();
			transfor_cor_best = transfor_cor;
		}
		if (score_best < 0.01) break;
	}
	//std::cout << "rcalTrans time:  " << calTrans_time << endl;
	//std::cout << "icp time:  " << icp_time << endl;

	pcl::transformPointCloud(*cloud_vol, *cloud_icp_trans, transfor_cor_icp);

	std::vector<int> nn_indices(1);
	std::vector<float> nn_dists(1);
	pcl::KdTreeFLANN<pcl::PointXYZ> tree_;
	tree_.setInputCloud(cloud_tgt_vol);
	// For each point in the source dataset
	int nr = 0;
	double fitness_score = 0;

	double max_range = 0.1;
	for (size_t i = 0; i < cloud_icp_trans->points.size(); ++i)
	{
		// Find its nearest neighbor in the target
		tree_.nearestKSearch(cloud_icp_trans->points[i], 1, nn_indices, nn_dists);

		// Deal with occlusions (incomplete targets)
		if (nn_dists[0] <= max_range)
		{
			//Add to the fitness score
			fitness_score += nn_dists[0];
			nr++;
		}

	}

	if (nr > 0)
		overlap = (double(nr) / cloud_tgt_vol->size());

	rough_output = transfor_cor_best;
	fine_output = transfor_cor_icp;
	return (score_best);

}