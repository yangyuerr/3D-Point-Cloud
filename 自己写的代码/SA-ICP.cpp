/*****************************************************************************
*                                SA-ICP.cpp
*作者毕业论文提出的算法
*
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/
double SA_ICP(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tgt, pcl::PointCloud<pcl::PointXYZ>::Ptr result, std::unordered_set<int> cyd_first,
	std::unordered_set<int> cyd_second, Eigen::Matrix4f transform_m, std::vector<std::vector<float>> transform_T)
{
	std::string mode_file_name = "database\\data.txt";
	std::ofstream mode_file(mode_file_name.c_str(), ios::out | ios::trunc);
	int iter = 100;
	//建立K-D tree
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree<pcl::PointXYZ>());
	tree_src->setInputCloud(cloud_tgt);
	//计算法向量
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal>ne_src;
	ne_src.setInputCloud(cloud_tgt);
	ne_src.setSearchMethod(tree_src);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_tgt_normals(new pcl::PointCloud<pcl::Normal>);
	ne_src.setRadiusSearch(Normal_dis);
	ne_src.compute(*cloud_tgt_normals);
	//首先建立最近点

	int n = 1000;
	float c_error = INT_MAX;
	Eigen::Matrix4f transform_total;
	//transform_total << 1, 0, 0, 0,
	//	0, 1, 0, 0,
	//	0, 0, 1, 0,
	//	0, 0, 0, 1;
	transform_total = transform_m;

	while (iter--)
	{
		//cout << iter << endl;
		pcl::PointCloud<pcl::PointXYZ>::Ptr nearest_points(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr source_points(new pcl::PointCloud<pcl::PointXYZ>);
		int n_sample = 1000;
		pcl::PointCloud<pcl::PointXYZ>::Ptr sample_cloud_src(new pcl::PointCloud<pcl::PointXYZ>);


		for (int i = 0; i < n_sample; ++i) {
			int k = 0;
			while (1)
			{
				k++;
				//	if (k >= 10000) break;
				int index = My_rand() % cloud_src->size();
				if (cyd_first.find(index) != cyd_first.end() || cyd_second.find(index) != cyd_second.end())
				{
					sample_cloud_src->push_back(cloud_src->at(index));
					break;
				}

			}
		}

		int num = 0;
		if (0)
		{
			for (int i = 0; i < cloud_src->size(); ++i)
			{
				std::vector<int> four_indices;
				std::vector<float> four_indices_dis;
				tree_src->nearestKSearch(*cloud_src, i, 1, four_indices, four_indices_dis);
				//if (sqrt(pow((*cloud_src)[i].x - cloud_tgt->at(four_indices[0]).x, 2) + pow((*sample_cloud_src)[i].y - cloud_tgt->at(four_indices[0]).y, 2) +
				//	pow((*sample_cloud_src)[i].z - cloud_tgt->at(four_indices[0]).z, 2))>2)
				//	continue;
				source_points->push_back(cloud_src->at(i));
				nearest_points->push_back(cloud_tgt->at(four_indices[0]));


			}
		}
		else
		{
			for (int i = 0; i < sample_cloud_src->size(); ++i)
			{
				std::vector<int> four_indices;
				std::vector<float> four_indices_dis;
				tree_src->nearestKSearch(*sample_cloud_src, i, 4, four_indices, four_indices_dis);
				pcl::PointCloud<pcl::Normal>::Ptr four_points_normals(new pcl::PointCloud<pcl::Normal>);
				pcl::PointCloud<pcl::PointXYZ>::Ptr four_point(new pcl::PointCloud<pcl::PointXYZ>);
				int _flag = 0;
				pcl::PointCloud<pcl::PointXYZ>::Ptr zzzz(new pcl::PointCloud<pcl::PointXYZ>);

				for (int j = 0; j < four_indices.size(); ++j)
				{
					//if (sqrt(pow((*cloud_src)[i].x - cloud_tgt->at(four_indices[j]).x, 2) + pow((*sample_cloud_src)[i].y - cloud_tgt->at(four_indices[j]).y, 2) +
					//	pow((*sample_cloud_src)[i].z - cloud_tgt->at(four_indices[j]).z, 2))>1)
					//if (abs(four_indices_dis[i])>0.015)
					//{
					//	num++;
					//	_flag = 1;
					//	break;

					//}
					//	else
					{
						//cout << abs(four_indices_dis[i]) << endl;
						four_points_normals->push_back(cloud_tgt_normals->at(four_indices[j]));
						four_point->push_back(cloud_tgt->at(four_indices[j]));
					}
				}
				//zzzz->push_back(sample_cloud_src->at(i));
				//visualize_pcd(zzzz, cloud_tgt,four_point);
				//if (_flag == 1) continue;


				////判断四个点是否属于一个圆柱
				//int _flag = 0;
				//if (cyd_first.find(four_indices[0]) != cyd_first.end() && cyd_first.find(four_indices[1]) != cyd_first.end() &&
				//	cyd_first.find(four_indices[2]) != cyd_first.end() && cyd_first.find(four_indices[3]) != cyd_first.end())
				//	_flag = 1;

				//if (cyd_second.find(four_indices[0]) != cyd_second.end() && cyd_second.find(four_indices[1]) != cyd_second.end() &&
				//	cyd_second.find(four_indices[2]) != cyd_second.end() && cyd_second.find(four_indices[3]) != cyd_second.end())
				//	_flag = 1;

				//if (_flag = 0) continue;



				pcl::ModelCoefficients::Ptr coeffients_cylinder(new pcl::ModelCoefficients);
				pcl::PointIndices::Ptr inliers_cylinder(new pcl::PointIndices);
				pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal> seg;
				seg.setOptimizeCoefficients(true);
				seg.setModelType(pcl::SACMODEL_CYLINDER);
				seg.setMethodType(pcl::SAC_RANSAC);
				////seg.setNormalDistanceWeight(0.1);
				//seg.setNormalDistanceWeight(0.2);		// for coarse data
				seg.setMaxIterations(1);
				//seg.setDistanceThreshold(0.3);
				//seg.setRadiusLimits(0, 6);				// radius is within 8 centimeters
				seg.setInputCloud(four_point);
				seg.setInputNormals(four_points_normals);

				seg.segment(*inliers_cylinder, *coeffients_cylinder);
				if (coeffients_cylinder->values.size() != 0)//能拟合一个圆柱
				{
					float x0 = coeffients_cylinder->values[0];
					float y0 = coeffients_cylinder->values[1];
					float z0 = coeffients_cylinder->values[2];
					float l = coeffients_cylinder->values[3];
					float m = coeffients_cylinder->values[4];
					float n = coeffients_cylinder->values[5];
					float r0 = coeffients_cylinder->values[6];
					double u = (sample_cloud_src->at(i).x - x0)*l + (sample_cloud_src->at(i).y - y0)*m + (sample_cloud_src->at(i).z - z0)*n;
					pcl::PointXYZ P_tmp;
					P_tmp.x = x0 + u*l;
					P_tmp.y = y0 + u*m;
					P_tmp.z = z0 + u*n;
					double distance_ = (sample_cloud_src->at(i).x - x0)*(sample_cloud_src->at(i).x - x0) +
						(sample_cloud_src->at(i).y - y0)*(sample_cloud_src->at(i).y - y0) +
						(sample_cloud_src->at(i).z - z0)*(sample_cloud_src->at(i).z - z0) - pow((l*(sample_cloud_src->at(i).x - x0)
						+ m*(sample_cloud_src->at(i).y - y0) + n*(sample_cloud_src->at(i).z - z0)), 2);
					distance_ = sqrt(abs(distance_));

					//cout << abs(distance_ - r0) << endl;
					if (abs(distance_ - r0) >10 || distance_ == 0) continue;
					pcl::PointXYZ nearnestpoint;
					nearnestpoint.x = P_tmp.x + (sample_cloud_src->at(i).x - P_tmp.x)*r0 / distance_;
					nearnestpoint.y = P_tmp.y + (sample_cloud_src->at(i).y - P_tmp.y)*r0 / distance_;
					nearnestpoint.z = P_tmp.z + (sample_cloud_src->at(i).z - P_tmp.z)*r0 / distance_;

					//if (sqrt(pow(nearnestpoint.x - cloud_tgt->at(four_indices[0]).x, 2) + pow(nearnestpoint.y - cloud_tgt->at(four_indices[0]).y, 2) +
					//	pow(nearnestpoint.z - cloud_tgt->at(four_indices[0]).z, 2)) > 0.1)
					//	continue;

					source_points->push_back(sample_cloud_src->at(i));
					nearest_points->push_back(nearnestpoint);
					//zzzz->push_back(nearnestpoint);
					//visualize_pcd(zzzz, cloud_tgt, four_point);
				}
			}
		}
		//cout << cloud_src->size() << endl;
		//cout << source_points->size() << endl;
		//计算中心点,为建立方程做准备
		//double m_x = 0, m_y = 0, m_z = 0;
		//double m_x_n = 0, m_y_n = 0, m_z_n = 0;
		//for (int i = 0; i < source_points->size(); ++i)
		//{
		//	m_x += source_points->at(i).x;
		//	m_y += source_points->at(i).y;
		//	m_z += source_points->at(i).z;
		//	m_x_n += nearest_points->at(i).x;
		//	m_y_n += nearest_points->at(i).y;
		//	m_z_n += nearest_points->at(i).z;
		//}
		//m_x /= source_points->size();
		//m_y /= source_points->size();
		//m_z /= source_points->size();
		//m_x_n /= source_points->size();
		//m_y_n /= source_points->size();
		//m_z_n /= source_points->size();

		//double A[3][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		//for (size_t i = 0; i < source_points->size(); ++i)//points.size()
		//{
		//	pcl::PointXYZ p_p = source_points->at(i);
		//	pcl::PointXYZ p_n = nearest_points->at(i);
		//	A[0][0] += (p_n.x - m_x_n)*(p_p.x - m_x);
		//	A[0][1] += (p_n.x - m_x_n)*(p_p.y - m_y);
		//	A[0][2] += (p_n.x - m_x_n)*(p_p.z - m_z);
		//	A[1][0] += (p_n.y - m_y_n)*(p_p.x - m_x);
		//	A[1][1] += (p_n.y - m_y_n)*(p_p.y - m_y);
		//	A[1][2] += (p_n.y - m_y_n)*(p_p.z - m_z);
		//	A[2][0] += (p_n.z - m_z_n)*(p_p.x - m_x);
		//	A[2][1] += (p_n.z - m_z_n)*(p_p.y - m_y);
		//	A[2][2] += (p_n.z - m_z_n)*(p_p.z - m_z);

		//}

		//double U[3][3];
		//double S[3];
		//double V[3][3];
		//int k = 3;
		//svd(A, k, U, S, V);
		////-----------------------计算R、T-----------------------//
		//double R[3][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		//double T[3] = { 0, 0, 0 };

		//for (int i = 0; i < 3; i++) {
		//	for (int j = 0; j < 3; j++) {
		//		R[i][j] = U[i][0] * V[0][j] + U[i][1] * V[1][j] + U[i][2] * V[2][j];
		//	}
		//}
		//T[0] = m_x_n - (R[0][0] * m_x + R[0][1] * m_y + R[0][2] * m_z);
		//T[1] = m_y_n - (R[1][0] * m_x + R[1][1] * m_y + R[1][2] * m_z);
		//T[2] = m_z_n - (R[2][0] * m_x + R[2][1] * m_y + R[2][2] * m_z);

		Eigen::Matrix4f transform_;
		//for (int i = 0; i < 3; i++) {
		//	for (int j = 0; j < 3; j++) {
		//		transform_(i, j) = R[i][j];
		//	}
		//}
		//transform_(0, 3) = T[0]; transform_(1, 3) = T[1]; transform_(2, 3) = T[2];
		//transform_(3, 0) = 0; transform_(3, 1) = 0; transform_(3, 2) = 0; transform_(3, 3) = 1;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tmp(new pcl::PointCloud<pcl::PointXYZ>);
		//pcl::transformPointCloud(*cloud_src, *cloud_tmp, transform_);
		//cloud_src = cloud_tmp;

		pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
		icp.setInputSource(source_points);
		icp.setInputTarget(nearest_points);
		// 最大迭代次数
		icp.setMaximumIterations(1);
		icp.align(*cloud_tmp);

		transform_ = icp.getFinalTransformation();
		pcl::transformPointCloud(*cloud_src, *cloud_tmp, transform_);
		cloud_src = cloud_tmp;
		transform_total = transform_* transform_total;
		std::vector<float> res;
		std::vector<std::vector<float>> transform_total_v = std::vector<std::vector<float>>(4, std::vector<float>(4, 0));;
		std::vector<std::vector<float>> transform_vT = std::vector<std::vector<float>>(4, std::vector<float>(4, 0));;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				//transform_total_v[i][j] = 0;
				transform_total_v[i][j] = transform_total(i, j);
				transform_vT[i][j] = transform_T[i][j];
			}
		}
		//transform_total_v[0][0] = 1;
		//transform_total_v[1][1] = 1;
		//transform_total_v[2][2] = 1;
		//transform_total_v[3][3] = 1;
		// transform_vT[0][3] -= 0.852123;
		// transform_vT[1][3] -= 0.019;
		// transform_vT[2][3] += 0.03525;
		//cout << transform_total_v[0][3] << "  " << transform_total_v[1][3] << "  " << transform_total_v[2][3] << endl;

		error_func(transform_total_v, transform_vT, res);
		//cout << "error_angle : " << res[0] << endl;
		//cout << "error_dis : " << res[1] << endl;
		cout << res[0] << "	" << res[1] << "	";

		std::vector<std::vector<float>> transform_E = std::vector<std::vector<float>>(4, std::vector<float>(4, 0));;
		std::vector<std::vector<float>> transform_vv = std::vector<std::vector<float>>(4, std::vector<float>(4, 0));;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				//transform_total_v[i][j] = 0;
				transform_E[i][j] = 0;
				transform_vv[i][j] = transform_(i, j);
			}
		}

		transform_E[0][0] = 1;
		transform_E[1][1] = 1;
		transform_E[2][2] = 1;
		transform_E[3][3] = 1;
		error_func(transform_vv, transform_E, res);

		
		//cout << "error_angle - iter : " << res[0] << endl;
		if (res[0] > 0 && res[0] < 180)
			float feafef = res[0];
		else
		{
			cout << endl;
			break;
		}
			
		//cout << "error_dis - iter: " << res[1] << endl;
		//计算残差
		double _error = 0;
		//for (int i = 0; i < source_points->size(); ++i)
		//{
		//	double tmp_diss = sqrt(fabs(pow(source_points->at(i).x - nearest_points->at(i).x, 2) + pow(source_points->at(i).y - nearest_points->at(i).y, 2) + pow(source_points->at(i).z - nearest_points->at(i).z, 2)));
		//	if (tmp_diss > 0)
		//	{
		//		_error += tmp_diss;
		//		//cout << "tmp_dis : " << tmp_diss << endl;
		//		//cout << "_error : " << _error << endl;
		//	}
		//}

		_error /= res[0];
		//if (_error == c_error) break;
		//cout << "error :" << '\t';
		//cout << _error << endl;
		//cout << c_error << endl;
		c_error = _error;
		//if (c_error == 0) { system("pause"); }
		//pcl::PointCloud<pcl::PointXYZ>::Ptr zzzzz(new pcl::PointCloud<pcl::PointXYZ>);
		//变换之前的采样点、目标点与变换之后的点
		//visualize_pcd(source_points, cloud_tgt, cloud_src);
		//变换之前的采样点最近点与目标点
		//visualize_pcd(source_points, cloud_tgt, nearest_points);
		//变换之前的采样点、目标点与变换之后的点
		//	visualize_pcd(source_points, cloud_tgt, cloud_src);

	}
	return 0;
}