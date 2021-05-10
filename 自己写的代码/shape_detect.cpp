/*****************************************************************************
*                                shape_detect.cpp
*作者自己写的利用RANSAC估计平面与圆柱
*
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/

#include "Pcl_Eigen.h"
#include "algorithm.h"
#include <unordered_set>
#include "svd.h"
/*
*产生0-300万的随机数
*/
int My_rand(){
	int a = rand() % 300;
	int b = rand() % 10000;
	return a * 10000 + b;
}

/*
*点云平面的检测
*输入需要检测的点云 输出检测出来的平面点与去点平面之后的点
*/

void detect_plain(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src, pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals,
	std::vector<int>& index_plain, std::vector<int>& index_no_plain)
{
    for (int index = 0; index < (*cloud_src).size(); index++)    { index_no_plain.push_back(index); }
	int m1 = 100;//检测m次疑似平面
	srand(time(NULL));
	while (m1)
	{
		m1--;
		int pointA_index = index_no_plain[My_rand() % index_no_plain.size()];
		pcl::PointXYZ pointA = (*cloud_src)[pointA_index];
		int m2 = 100;
		while (m2)
		{
			m2--;
			int pointB_index = index_no_plain[My_rand() % index_no_plain.size()];
			pcl::PointXYZ pointB = (*cloud_src)[pointB_index];
			double normals_tmp = sqrt((pointB.x - pointA.x)*(pointB.x - pointA.x) + (pointB.y - pointA.y)*(pointB.y - pointA.y)
				+ (pointB.z - pointA.z)*(pointB.z - pointA.z));
			double normalsA_AB = (((pointB.x - pointA.x)*(*cloud_src_normals)[pointA_index].normal_x)
				+ ((pointB.y - pointA.y)*(*cloud_src_normals)[pointA_index].normal_y)
				+ ((pointB.z - pointA.z)*(*cloud_src_normals)[pointA_index].normal_z)) / normals_tmp;
			double normalsB_AB = (((pointB.x - pointA.x)*(*cloud_src_normals)[pointB_index].normal_x)
				+ ((pointB.y - pointA.y)*(*cloud_src_normals)[pointB_index].normal_y)
				+ ((pointB.z - pointA.z)*(*cloud_src_normals)[pointB_index].normal_z)) / normals_tmp;
			if (abs(normalsA_AB)>cos_angle_90err_plain || abs(normalsB_AB)>cos_angle_90err_plain)
			{
				continue;
			}
			else
			{
				int m3 = 100;
				while (m3){
					m3--;
					int pointC_index = index_no_plain[My_rand() % index_no_plain.size()];
					pcl::PointXYZ pointC = (*cloud_src)[pointC_index];
					normals_tmp = sqrt((pointC.x - pointA.x)*(pointC.x - pointA.x) + (pointC.y - pointA.y)*(pointC.y - pointA.y)
						+ (pointC.z - pointA.z)*(pointC.z - pointA.z));
					double normalsA_AC = (((pointC.x - pointA.x)*(*cloud_src_normals)[pointA_index].normal_x)
						+ ((pointC.y - pointA.y)*(*cloud_src_normals)[pointA_index].normal_y)
						+ ((pointC.z - pointA.z)*(*cloud_src_normals)[pointA_index].normal_z)) / normals_tmp;
					double normalsC_AC = (((pointC.x - pointA.x)*(*cloud_src_normals)[pointC_index].normal_x)
						+ ((pointC.y - pointA.y)*(*cloud_src_normals)[pointC_index].normal_y)
						+ ((pointC.z - pointA.z)*(*cloud_src_normals)[pointC_index].normal_z)) / normals_tmp;

					normals_tmp = sqrt((pointC.x - pointB.x)*(pointC.x - pointB.x) + (pointC.y - pointB.y)*(pointC.y - pointB.y)
						+ (pointC.z - pointB.z)*(pointC.z - pointB.z));
					double normalsB_BC = (((pointC.x - pointB.x)*(*cloud_src_normals)[pointB_index].normal_x)
						+ ((pointC.y - pointB.y)*(*cloud_src_normals)[pointB_index].normal_y)
						+ ((pointC.z - pointB.z)*(*cloud_src_normals)[pointB_index].normal_z)) / normals_tmp;
					double normalsC_BC = (((pointC.x - pointB.x)*(*cloud_src_normals)[pointC_index].normal_x)
						+ ((pointC.y - pointB.y)*(*cloud_src_normals)[pointC_index].normal_y)
						+ ((pointC.z - pointB.z)*(*cloud_src_normals)[pointC_index].normal_z)) / normals_tmp;

					//计算平面
					double A = (pointC.y - pointA.y)*(pointC.z - pointA.z) - (pointB.z - pointA.z)*(pointC.y - pointA.y);
					double B = (pointC.x - pointA.x)*(pointB.z - pointA.z) - (pointB.x - pointA.x)*(pointC.z - pointA.z);
					double C = (pointB.x - pointA.x)*(pointC.y - pointA.y) - (pointC.x - pointA.x)*(pointB.y - pointA.y);
					double D = -(A*pointA.x + B*pointA.y + C*pointA.z);

					double temp_angle_dis = sqrt(A*A + B*B + C*C);
					double cosA = (A*((*cloud_src_normals)[pointA_index]).normal_x + B*((*cloud_src_normals)[pointA_index]).normal_y
						+ C*((*cloud_src_normals)[pointA_index]).normal_z) / temp_angle_dis;
					double cosB = (A*((*cloud_src_normals)[pointB_index]).normal_x + B*((*cloud_src_normals)[pointB_index]).normal_y
						+ C*((*cloud_src_normals)[pointB_index]).normal_z) / temp_angle_dis;
					double cosC = (A*((*cloud_src_normals)[pointC_index]).normal_x + B*((*cloud_src_normals)[pointC_index]).normal_y
						+ C*((*cloud_src_normals)[pointC_index]).normal_z) / temp_angle_dis;

					if (abs(normalsA_AC) > cos_angle_90err_plain || abs(normalsC_AC) > cos_angle_90err_plain || abs(normalsB_BC) > cos_angle_90err_plain || abs(normalsC_BC) > cos_angle_90err_plain
						|| (1 - abs(cosA))>cos_angle_0err_plain || (1 - abs(cosB))>cos_angle_0err_plain || (1 - abs(cosC))>cos_angle_0err_plain)
					{
						continue;
					}
					else
					{
						//计算内点
						double temp_dis = sqrt(A*A+B*B+C*C+D*D);
						std::vector<int> index_plain_tmp;
						std::vector<int> index_no_plain_tmp;
						for (int i = 0; i < index_no_plain.size(); i++){
							double distance_ = (((*cloud_src)[index_no_plain[i]].x)*A + ((*cloud_src)[index_no_plain[i]].y)*B +
								((*cloud_src)[index_no_plain[i]].z)*C + D) / temp_dis;
							distance_ = abs(distance_);

							double cos_angle = (A*(*cloud_src_normals)[index_no_plain[i]].normal_x +
								B*(*cloud_src_normals)[index_no_plain[i]].normal_y +
								C*(*cloud_src_normals)[index_no_plain[i]].normal_z) / temp_angle_dis;

							double sin_angle = sqrt(1 - cos_angle*cos_angle);
							if (distance_ < max_D&&(1-abs(cos_angle)<cos_angle_0err_plain)){
								index_plain_tmp.push_back(index_no_plain[i]);
							}
							else{
								index_no_plain_tmp.push_back(index_no_plain[i]);
							}
						}

						index_no_plain = index_no_plain_tmp;
						for (int i = 0; i < index_plain_tmp.size(); i++){ index_plain.push_back(index_plain_tmp[i]); }
						m3 = 0;
						m2 = 0;
					}

				}

			}

		}

	}
}

/*
*此函数用来判断两个点是否可能是一个圆柱上的点
*/
bool TwopointsCyd(PointT& A, PointT& B, pcl::Normal& N_a, pcl::Normal& N_b, float& radius, PointT *point2cenA, PointT *point2cenB, int flag)
{
	int m = 4;
	while (m)
	{
		m--;
		if (m == 3)
		{
			if (flag == 0)
			{
				point2cenA->x = A.x + radius*N_a.normal_x;
				point2cenA->y = A.y + radius*N_a.normal_y;
				point2cenA->z = A.z + radius*N_a.normal_z;
				point2cenB->x = B.x + radius*N_b.normal_x;
				point2cenB->y = B.y + radius*N_b.normal_y;
				point2cenB->z = B.z + radius*N_b.normal_z;
			}
			else if (flag == 1)
			{
				point2cenB->x = B.x + radius*N_b.normal_x;
				point2cenB->y = B.y + radius*N_b.normal_y;
				point2cenB->z = B.z + radius*N_b.normal_z;
			}

		}
		else if (m == 2)
		{
			if (flag == 0)
			{
				point2cenA->x = A.x - radius*N_a.normal_x;
				point2cenA->y = A.y - radius*N_a.normal_y;
				point2cenA->z = A.z - radius*N_a.normal_z;
				point2cenB->x = B.x - radius*N_b.normal_x;
				point2cenB->y = B.y - radius*N_b.normal_y;
				point2cenB->z = B.z - radius*N_b.normal_z;
			}
			else if (flag == 1)
			{
				point2cenB->x = B.x - radius*N_b.normal_x;
				point2cenB->y = B.y - radius*N_b.normal_y;
				point2cenB->z = B.z - radius*N_b.normal_z;
			}

		}
		else if (m == 1)
		{
			if (flag == 0)
			{
				point2cenA->x = A.x - radius*N_a.normal_x;
				point2cenA->y = A.y - radius*N_a.normal_y;
				point2cenA->z = A.z - radius*N_a.normal_z;
				point2cenB->x = B.x + radius*N_b.normal_x;
				point2cenB->y = B.y + radius*N_b.normal_y;
				point2cenB->z = B.z + radius*N_b.normal_z;
			}
			else if (flag == 1)
			{
				point2cenB->x = B.x + radius*N_b.normal_x;
				point2cenB->y = B.y + radius*N_b.normal_y;
				point2cenB->z = B.z + radius*N_b.normal_z;
			}

		}
		else
		{
			if (flag == 0)
			{
				point2cenA->x = A.x + radius*N_a.normal_x;
				point2cenA->y = A.y + radius*N_a.normal_y;
				point2cenA->z = A.z + radius*N_a.normal_z;
				point2cenB->x = B.x - radius*N_b.normal_x;
				point2cenB->y = B.y - radius*N_b.normal_y;
				point2cenB->z = B.z - radius*N_b.normal_z;
			}
			else if (flag == 1)
			{
				point2cenB->x = B.x - radius*N_b.normal_x;
				point2cenB->y = B.y - radius*N_b.normal_y;
				point2cenB->z = B.z - radius*N_b.normal_z;
			}
		}

		float Normal_tmp = sqrt((point2cenA->x - point2cenB->x)*(point2cenA->x - point2cenB->x) +
			(point2cenA->y - point2cenB->y)*(point2cenA->y - point2cenB->y)
			+ (point2cenA->z - point2cenB->z)*(point2cenA->z - point2cenB->z));
		if (Normal_tmp < 0.4) return true;
		float normalA_AB = ((point2cenA->x - point2cenB->x)*N_a.normal_x +
			(point2cenA->y - point2cenB->y)*N_a.normal_y +
			(point2cenA->z - point2cenB->z)*N_a.normal_z) / Normal_tmp;
		float normalB_AB = ((point2cenA->x - point2cenB->x)*N_b.normal_x +
			(point2cenA->y - point2cenB->y)*N_b.normal_y +
			(point2cenA->z - point2cenB->z)*N_b.normal_z) / Normal_tmp;
		if (fabs(normalA_AB) <= cos_angle_90err_cyd&&fabs(normalB_AB) <= cos_angle_90err_cyd)
		{
			return true;
		}
	}
	return false;
}

/*
*多次迭代圆柱拟合圆柱
*输入原始点云 初始圆柱内点 可用候选点（有内点之后 圆柱候选点） 初始参数与最终参数
*圆柱轴线的长度与中心点的位置，默认输入 iter 迭代次数 默认为4
*输出是否可以拟合一个圆柱
*
*/
bool cylinder_fit_iter(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src, std::vector<int>& index_cyd, std::vector<int>& index,
	cylinder_coeffient_t *initial_params, cylinder_coeffient_t *final_params, std::vector<float>& axis_length, int iter = 4)
{
	
	int points = index_cyd.size();
	while (iter--){
		std::vector<int> index_cyd_tmp;
		std::vector<int> index_no_cyd_tmp;
		std::vector<float> data_dtr;
		std::vector<int> index_cyd_cyd;
		float min_point = INT_MAX;
		int index_min_point=0;
		float max_point = INT_MIN;
		int index_max_point=0;
		data_dtr = cylinder_fit(cloud_src, index_cyd, initial_params, final_params);

		//做第一层筛选，去掉一些有问题的点
		for (int i = 0; i < index_cyd.size(); i++){
			double x = (*cloud_src)[index_cyd[i]].x;
			double y = (*cloud_src)[index_cyd[i]].y;
			double z = (*cloud_src)[index_cyd[i]].z;
			double distance_ = (x - final_params->x)*(x - final_params->x) +
				(y - final_params->y)*(y - final_params->y) +
				(z - final_params->z)*(z - final_params->z) - pow((final_params->a*(x - final_params->x)
				+ final_params->b*(y - final_params->y) + final_params->c*(z - final_params->z)), 2);
			distance_ = sqrt(abs(distance_));
			distance_ = abs(distance_ - final_params->r);
			if (distance_ < data_dtr[0])
			{
				index_cyd_cyd.push_back(index_cyd[i]);
			}
		}
		data_dtr = cylinder_fit(cloud_src, index_cyd_cyd, final_params, final_params);
		for (int i = 0; i < index.size(); i++){
			double x = (*cloud_src)[index[i]].x;
			double y = (*cloud_src)[index[i]].y;
			double z = (*cloud_src)[index[i]].z;
			double distance_ = (x - final_params->x)*(x - final_params->x) +
				(y - final_params->y)*(y - final_params->y) +
				(z - final_params->z)*(z - final_params->z) - pow((final_params->a*(x - final_params->x)
				+ final_params->b*(y - final_params->y) + final_params->c*(z - final_params->z)), 2);
			distance_ = sqrt(abs(distance_));
			distance_ = abs(distance_ - final_params->r);

			if (distance_ < 0.05){
				float lgh = ((x - final_params->x)*(final_params->a) + (y - final_params->y)*final_params->b +
					(z - final_params->z)*(final_params->c)) / (sqrt((final_params->a)*(final_params->a) +
					(final_params->b)*(final_params->b) + (final_params->c)*(final_params->c)));

				if (min_point > lgh) { min_point = lgh; index_min_point = index[i]; }
				if (max_point < lgh) { max_point = lgh; index_max_point = index[i]; }
				//min_point = std::min(min_point, lgh);
				//max_point = std::max(max_point, lgh);
				index_cyd_tmp.push_back(index[i]);
			}
			else{
				index_no_cyd_tmp.push_back(index[i]);
			}
		}

		index_cyd = index_cyd_tmp;
		initial_params = final_params;
		axis_length = { max_point, min_point,float(index_max_point),float(index_min_point) };

		if (index_cyd_tmp.size() == 0) return false;
		if (points / index_cyd_tmp.size() >= 2)  return false;
		else if (points == index_cyd_tmp.size()) break;
		points = index_cyd_tmp.size();
	}
	return true;
}




/*
*最小二乘法拟合圆柱
*
*/
std::vector<float> cylinder_fit(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in, std::vector<int>& index, cylinder_coeffient_t *initial_params, cylinder_coeffient_t *final_params)
{
	using namespace Eigen;
	MatrixXf params(6, 1);
	params << 0, 0, 0, 0, 0, 0;
	MatrixXf X(6, 1);
	X << initial_params->x, initial_params->y, initial_params->z,
		initial_params->a, initial_params->b, initial_params->c;
	int times = 0;
	int size = index.size();
	MatrixXf B(size, 6), BT(6, size), L(size, 1);
	while ((fabs(X(0, 0)) + fabs(X(1, 0)) + fabs(X(2, 0)) > 0.1) || times == 0)
	{
		if (times >= 10) break;
		params += X;
		float x0 = params(0, 0);
		float y0 = params(1, 0);
		float z0 = params(2, 0);
		float a = params(3, 0);
		float b = params(4, 0);
		float c = params(5, 0);
		float r = initial_params->r;
		float f0;
		for (int i = 0; i < index.size(); i++)
		{

			float x = (*cloud_in)[index[i]].x;
			float y = (*cloud_in)[index[i]].y;
			float z = (*cloud_in)[index[i]].z;

			f0 = pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2) - pow((a*(x - x0) + b * (y - y0) + c * (z - z0)), 2) - r * r;
			L(i, 0) = -f0;
			B(i, 0) = (a * (a*(x - x0) + b * (y - y0) + c * (z - z0)) - (x - x0)) * 2;
			B(i, 1) = (b * (a*(x - x0) + b * (y - y0) + c * (z - z0)) - (y - y0)) * 2;
			B(i, 2) = (c * (a*(x - x0) + b * (y - y0) + c * (z - z0)) - (z - z0)) * 2;
			B(i, 3) = -(a*(x - x0) + b * (y - y0) + c * (z - z0))*(x - x0) * 2;
			B(i, 4) = -(a*(x - x0) + b * (y - y0) + c * (z - z0))*(y - y0) * 2;
			B(i, 5) = -(a*(x - x0) + b * (y - y0) + c* (z - z0))*(z - z0) * 2;

		}
		X = B.colPivHouseholderQr().solve(L);
		times++;
	}
	final_params->x = params(0, 0);
	final_params->y = params(1, 0);
	final_params->z = params(2, 0);
	final_params->a = params(3, 0);
	final_params->b = params(4, 0);
	final_params->c = params(5, 0);
	final_params->r = initial_params->r;
	std::vector<float> res;
	float mean = 0;
	float var = 0;
	for (int i = 0; i < index.size(); i++){
		double x = (*cloud_in)[index[i]].x;
		double y = (*cloud_in)[index[i]].y;
		double z = (*cloud_in)[index[i]].z;
		double distance_ = (x - final_params->x)*(x - final_params->x) +
			(y - final_params->y)*(y - final_params->y) +
			(z - final_params->z)*(z - final_params->z) - pow((final_params->a*(x - final_params->x)
			+ final_params->b*(y - final_params->y) + final_params->c*(z - final_params->z)), 2);
		distance_ = sqrt(abs(distance_));
		distance_ = abs(distance_ - final_params->r);

		mean += distance_;
		var += distance_*distance_;
	}
	res.push_back(mean / index.size());
	res.push_back(var / index.size());
	return res;
}

/*
*计算圆柱点云的覆盖率

*/
double cover_rate_cyd(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src, std::vector<int>& index, cylinder_coeffient_t * final_params, std::vector<float> axis_length)
{
	//计算经过最小点并垂直于圆柱轴线的平面，然后计算所有的点到这个平面的距离，划分区间，然后计算每个区间的点数，
	//当小于某个值就认为该区间不计算在内，之后计算一个覆盖率 总点数（需要将可行区间内的点算入）*分辨率/pi * r2*l
	int H_size = int((axis_length[0] - axis_length[1]) / ratio) + 2;
	
	std::vector<int> H = std::vector<int>(H_size, 0);
	float A = final_params->a;
	float B = final_params->b;
	float C = final_params->c;
	float r = final_params->r;
	int index_min_point = int(axis_length[3]);
	float D = -(A*(*cloud_src)[index_min_point].x + B*(*cloud_src)[index_min_point].y + C*(*cloud_src)[index_min_point].z);
	for (int i = 0; i < index.size(); ++i)
	{
		pcl::PointXYZ p_tmp = (*cloud_src)[index[i]];
		float dis_ = (A*p_tmp.x + B*p_tmp.y + C*p_tmp.z + D) / sqrt(A*A+B*B+C*C);
		H[int(dis_ / ratio)] += 1;

	}
	int L = 0;
	int points_num = 0;
	for (int i = 0; i < H_size; ++i)
	{
		if (H[i] > 10)
		{
			L += 1;
			points_num += H[i];
		}
	}
	if (L == 0) return 0;
	return points_num*ratio / (M_PI*2*r*L);

}
/*
*圆柱检测
*输入需要检测的点云 输出检测出来的圆柱及其参数 以及圆柱轴线线段
*/
void detect_cyd(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src, std::vector<int>& index, pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals,
std::vector<int>& syd_candidate,std::vector<cylinder_coeffient_t*>& cyd_cof_candi,
std::vector<std::vector<float>>& axis_length_candi)
{
    srand(time(NULL));
	std::vector<int> Data_set = index;//数据集
    std::vector<float> radius = { float(2.05), float(1.325)};
	int m1 = 500;
	pcl::PointCloud<pcl::PointXYZ>::Ptr zzzzz(new pcl::PointCloud<pcl::PointXYZ>);
	while (m1){
		m1--;
		//cout << m1 << endl;
		int pointA_index = Data_set[ My_rand() % Data_set.size()];
		pcl::PointXYZ pointA = (*cloud_src)[pointA_index];
		//具体哪个候选半径
		int radius_candidate = 0;
		//用来记录大致的轴向
		std::vector<float> Normal_candidate;
		pcl::PointCloud<pcl::PointXYZ>::Ptr first(new pcl::PointCloud<pcl::PointXYZ>);
		
		int m2 = 100;
		while (m2){
			m2--;
			int pointB_index = Data_set[My_rand() % Data_set.size()];
			pcl::PointXYZ pointB = (*cloud_src)[pointB_index];
			float tmp_cos=(*cloud_src_normals)[pointB_index].normal_x*(*cloud_src_normals)[pointA_index].normal_x+
				(*cloud_src_normals)[pointB_index].normal_y*(*cloud_src_normals)[pointA_index].normal_y + 
				(*cloud_src_normals)[pointB_index].normal_z*(*cloud_src_normals)[pointA_index].normal_z;
			if ((1 - abs(tmp_cos)) <= cos_angle_0err_cyd_sample_normals) continue; //去掉平面和线上的点

			float tmp_dis = pow(pointB.x-pointA.x,2) + pow(pointA.y-pointB.y,2) + pow(pointA.z*pointB.z,2);//去掉远距离的点
			if (tmp_dis>=9) continue;


			PointT* point2cenA=new PointT;
			PointT* point2cenB = new PointT;
			bool flag = TwopointsCyd(pointA, pointB, (*cloud_src_normals)[pointA_index], (*cloud_src_normals)[pointB_index],
				radius[0], point2cenA, point2cenB, 0);
			if (flag == 1) radius_candidate = 1;//第一个半径
			else
			{
				flag = TwopointsCyd(pointA, pointB, (*cloud_src_normals)[pointA_index], (*cloud_src_normals)[pointB_index],
					radius[1], point2cenA, point2cenB, 0);
				if (flag == 1) radius_candidate = 2;//第二个半径
				else continue;
			}
			int m3 = 100;
			while (m3){
				m3--;
				int pointC_index = Data_set[My_rand() % Data_set.size()];
				pcl::PointXYZ pointC = (*cloud_src)[pointC_index];
				tmp_cos = (*cloud_src_normals)[pointB_index].normal_x*(*cloud_src_normals)[pointC_index].normal_x +
					(*cloud_src_normals)[pointB_index].normal_y*(*cloud_src_normals)[pointC_index].normal_y +
					(*cloud_src_normals)[pointB_index].normal_z*(*cloud_src_normals)[pointC_index].normal_z;
				if ((1 - abs(tmp_cos)) <= cos_angle_0err_cyd_sample_normals) continue; //去掉平面和线上的点

				tmp_cos = (*cloud_src_normals)[pointC_index].normal_x*(*cloud_src_normals)[pointA_index].normal_x +
					(*cloud_src_normals)[pointC_index].normal_y*(*cloud_src_normals)[pointA_index].normal_y +
					(*cloud_src_normals)[pointC_index].normal_z*(*cloud_src_normals)[pointA_index].normal_z;
				if ((1 - abs(tmp_cos)) <= cos_angle_0err_cyd_sample_normals) continue; //去掉平面和线上的点

				tmp_dis = pow(pointC.x-pointA.x,2) + pow(pointA.y-pointC.y,2) + pow(pointA.z-pointC.z,2);//去掉远距离的点
				if (tmp_dis >= 9) continue;
				tmp_dis = pow(pointB.x-pointC.x,2) + pow(pointC.y-pointB.y,2) + pow(pointC.z-pointB.z,2);//去掉远距离的点
				if (tmp_dis >= 9) continue;


				PointT* point2cenC = new PointT;

				flag = TwopointsCyd(pointA, pointC, (*cloud_src_normals)[pointA_index], (*cloud_src_normals)[pointC_index],
					radius[radius_candidate - 1], point2cenA, point2cenC, 1);
				if (flag)
				{
					flag = TwopointsCyd(pointB, pointC, (*cloud_src_normals)[pointB_index], (*cloud_src_normals)[pointC_index],
						radius[radius_candidate - 1], point2cenA, point2cenC, 2);

					if (flag){
						first->push_back(pointA);
						first->push_back(pointB);
						first->push_back(pointC);
						first->push_back(*point2cenA);
						first->push_back(*point2cenB);
						first->push_back(*point2cenC);
						//三个点、一个半径、大致轴向 可以确定一个大致的圆柱
						//从所有的点中做一个初步的筛选，得到第一次的内点
						//根据内点得到一个接近精确地圆柱方程
						//再从所有的点中得到在这个圆柱上的点
						std::vector<int> index_cyd;
						Normal_candidate.push_back(point2cenB->x + point2cenC->x - 2 * point2cenA->x);
						Normal_candidate.push_back(point2cenB->y + point2cenC->y - 2 * point2cenA->y);
						Normal_candidate.push_back(point2cenB->z + point2cenC->z - 2 * point2cenA->z);
						PointT* point2cenD = new PointT;
						PointT tmp_Point;
						tmp_Point.x = tmp_Point.y = tmp_Point.z = 0;
						int times = 0;
						for (int i = 0; i < Data_set.size(); i++)
						{
							pcl::PointXYZ pointD = (*cloud_src)[Data_set[i]];
							bool flag1 = TwopointsCyd(pointA, pointD, (*cloud_src_normals)[pointA_index], (*cloud_src_normals)[Data_set[i]],
								radius[radius_candidate - 1], point2cenA, point2cenD, 1);
							bool flag2 = TwopointsCyd(pointB, pointD, (*cloud_src_normals)[pointB_index], (*cloud_src_normals)[Data_set[i]],
								radius[radius_candidate - 1], point2cenB, point2cenD, 2);
							bool flag3 = TwopointsCyd(pointC, pointD, (*cloud_src_normals)[pointC_index], (*cloud_src_normals)[Data_set[i]],
								radius[radius_candidate - 1], point2cenC, point2cenD, 2);
							if (flag1&&flag2&&flag3){
								times++;
								tmp_Point.x += pointD.x;
								tmp_Point.y += pointD.y;
								tmp_Point.z += pointD.z;
								index_cyd.push_back(Data_set[i]);

							}
						}

						float tmp_normals = sqrt(Normal_candidate[0] * Normal_candidate[0] + Normal_candidate[1] * Normal_candidate[1] +
							Normal_candidate[2] * Normal_candidate[2]);
						cylinder_coeffient_t * initial_params =new cylinder_coeffient_t;
						initial_params->x = tmp_Point.x / times;
						initial_params->y = tmp_Point.y / times;
						initial_params->z = tmp_Point.z / times;
						initial_params->a = Normal_candidate[0] / tmp_normals;
						initial_params->b = Normal_candidate[1] / tmp_normals;
						initial_params->c = Normal_candidate[2] / tmp_normals;
						initial_params->r = radius[radius_candidate - 1];
						cylinder_coeffient_t * final_params = new cylinder_coeffient_t;
						std::vector<float> axis_length;

						if (cylinder_fit_iter(cloud_src, index_cyd, Data_set, initial_params, final_params, axis_length)&&axis_length.size()!=0)
						{
							//计算圆柱的点云的覆盖率，将点云覆盖率小于某个阈值的删掉
							float rate = cover_rate_cyd(cloud_src, index_cyd, final_params, axis_length);
							if (rate<0.3) continue;
							for (int i = 0; i < index_cyd.size(); ++i) zzzzz->push_back((*cloud_src)[index_cyd[i]]);

							syd_candidate.push_back(index_cyd.size());
							cyd_cof_candi.push_back(final_params);
							axis_length_candi.push_back(axis_length);
						}
						m3 = 0;
						m2 = 0;
					}
				}
			}
		}
	}
	pcl::PointCloud<pcl::PointXYZ>::Ptr aaa(new pcl::PointCloud<pcl::PointXYZ>);
	//visualize_pcd(aaa, zzzzz, cloud_src);

}

/*
*此函数用来求空间中两条线的距离与最近点
*输入为直线A与直线B上的两个点
*输出为距离与最近的两个点
*
*/
void line2line_dis(pcl::PointXYZ a1, pcl::PointXYZ a2, pcl::PointXYZ b1, pcl::PointXYZ b2, std::vector<double>& res)
{
	pcl::PointXYZ v1 = Pointminus(a1, b1), v2 = Pointminus(a2, b2);
	pcl::PointXYZ N1 = PointCross(v1, v2);
	pcl::PointXYZ ab = Pointminus(a1, a2);
	double ans = PointDot(N1, ab) / Pointlength(N1);
	pcl::PointXYZ p1 = a1, p2 = a2;
	pcl::PointXYZ d1 = Pointminus(b1, a1), d2 = Pointminus(b2, a2);
	pcl::PointXYZ ans1, ans2;
	double t1, t2;
	t1 = PointDot((PointCross(Pointminus(p2, p1), d2)), PointCross(d1, d2));
	t2 = PointDot((PointCross(Pointminus(p2, p1), d1)), PointCross(d1, d2));
	double dd = Pointlength((PointCross(d1, d2)));
	t1 /= dd * dd;
	t2 /= dd * dd;
	ans1 = Pointadd(a1, Pointmulti(Pointminus(b1, a1), t1));
	ans2 = Pointadd(a2, Pointmulti(Pointminus(b2, a2), t2));
	res.push_back(ans);
	res.push_back(ans1.x);
	res.push_back(ans1.y);
	res.push_back(ans1.z);
	res.push_back(ans2.x);
	res.push_back(ans2.y);
	res.push_back(ans2.z);
}

/*
*此函数用来两个轴向计算RT方程
*
*/
void twoaxis2trans(pcl::PointXYZ center, pcl::Normal a, pcl::Normal b, Eigen::Affine3f& transform)
{

	pcl::PointXYZ center_target;
	center_target.x = 0;
	center_target.y = 0;
	center_target.z = 0;
	pcl::Normal a_target;
	a_target.normal_x = -1;
	a_target.normal_y = 0;
	a_target.normal_z = 0;
	pcl::Normal b_target;
	b_target.normal_x = 0;
	b_target.normal_y = 1;
	b_target.normal_z = 0;

	float dis_a = sqrt(a.normal_x*a.normal_x + a.normal_y*a.normal_y + a.normal_z*a.normal_z);
	a.normal_x = a.normal_x / dis_a;
	a.normal_y = a.normal_y / dis_a;
	a.normal_z = a.normal_z / dis_a;

	float dis_b = sqrt(b.normal_x*b.normal_x + b.normal_y*b.normal_y + b.normal_z*b.normal_z);
	b.normal_x = b.normal_x / dis_b;
	b.normal_y = b.normal_y / dis_b;
	b.normal_z = b.normal_z / dis_b;
	Eigen::Affine3f transform_T = Eigen::Affine3f::Identity();
	//先平移，对齐中心点
	transform_T.translation() << center_target.x - center.x, center_target.y - center.y, center_target.z - center.z;


	//先还原（-1，0，0）朝向,T字的一竖
	//先绕Z轴转动到XZ平面,再绕Y轴旋转-1 0 0
	Eigen::Affine3f transform_R1 = Eigen::Affine3f::Identity();
	Eigen::Affine3f transform_R2 = Eigen::Affine3f::Identity();
	float cos_thetaz = a.normal_x / sqrt(a.normal_x*a.normal_x + a.normal_y*a.normal_y);
	float sin_thetaz;
	if (a.normal_y < 0)
		sin_thetaz = sqrt(1 - cos_thetaz * cos_thetaz);
	else
		sin_thetaz = -sqrt(1 - cos_thetaz * cos_thetaz);


	float cos_thetay = -sqrt(a.normal_x*a.normal_x + a.normal_y*a.normal_y);
	//cos_thetay = a.normal_x<0 ? cos_thetay : -cos_thetay;
	float sin_thetay;
	if (a.normal_z < 0)
		sin_thetay = sqrt(1 - cos_thetay * cos_thetay);
	else
		sin_thetay = -sqrt(1 - cos_thetay * cos_thetay);

	transform_R1(0, 0) = cos_thetaz;
	transform_R1(0, 1) = -sin_thetaz;
	transform_R1(0, 2) = 0;
	transform_R1(1, 0) = sin_thetaz;
	transform_R1(1, 1) = cos_thetaz;
	transform_R1(1, 2) = 0;
	transform_R1(2, 0) = 0;
	transform_R1(2, 1) = 0;
	transform_R1(2, 2) = 1;
	transform_R2(0, 0) = cos_thetay;
	transform_R2(0, 1) = 0;
	transform_R2(0, 2) = sin_thetay;
	transform_R2(1, 0) = 0;
	transform_R2(1, 1) = 1;
	transform_R2(1, 2) = 0;
	transform_R2(2, 0) = -sin_thetay;
	transform_R2(2, 1) = 0;
	transform_R2(2, 2) = cos_thetay;


	pcl::Normal b_after1;

	b_after1.normal_x = transform_R1(0, 0) * b.normal_x + transform_R1(0, 1) * b.normal_y + transform_R1(0, 2) * b.normal_z;
	b_after1.normal_y = transform_R1(1, 0) * b.normal_x + transform_R1(1, 1) * b.normal_y + transform_R1(1, 2) * b.normal_z;
	b_after1.normal_z = transform_R1(2, 0) * b.normal_x + transform_R1(2, 1) * b.normal_y + transform_R1(2, 2) * b.normal_z;
	//cout << b_after1.normal_x << '\t' << b_after1.normal_y << '\t' << b_after1.normal_z << endl;

	pcl::Normal b_after2;

	b_after2.normal_x = transform_R2(0, 0) * b_after1.normal_x + transform_R2(0, 1) * b_after1.normal_y + transform_R2(0, 2) * b_after1.normal_z;
	b_after2.normal_y = transform_R2(1, 0) * b_after1.normal_x + transform_R2(1, 1) * b_after1.normal_y + transform_R2(1, 2) * b_after1.normal_z;
	b_after2.normal_z = transform_R2(2, 0) * b_after1.normal_x + transform_R2(2, 1) * b_after1.normal_y + transform_R2(2, 2) * b_after1.normal_z;
	//cout << b_after2.normal_x << '\t' << b_after2.normal_y << '\t' << b_after2.normal_z << endl;

	//在还原（0，1，0）朝向,T字的一横
	Eigen::Affine3f transform_R3 = Eigen::Affine3f::Identity();

	float cos_thetax = b_after2.normal_y / sqrt(b_after2.normal_y*b_after2.normal_y + b_after2.normal_z*b_after2.normal_z);
	float sin_thetax;
	if (b_after2.normal_z > 0)
		sin_thetax = sqrt(1 - cos_thetax * cos_thetax);
	else
		sin_thetax = -sqrt(1 - cos_thetax * cos_thetax);

	transform_R3(1, 1) = cos_thetax;
	transform_R3(1, 2) = sin_thetax;
	transform_R3(2, 1) = -sin_thetax;
	transform_R3(2, 2) = cos_thetax;

	transform = transform_R3 *transform_R2* transform_R1 * transform_T;
}


/*
*根据检测的圆柱进行匹配
*
*/
void registration(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_src,pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_tgt,
std::vector<int>& syd_candidate,std::vector<cylinder_coeffient_t*>& cyd_cof_candi,
std::vector<std::vector<float>>& axis_length_candi, std::vector<std::vector<float>>& Transform,
pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_final, std::vector<cylinder_coeffient_t*>& cyd_vertical)
{

	//pcl::PointCloud<pcl::PointXYZ>::Ptr zzzzz(new pcl::PointCloud<pcl::PointXYZ>);
	//pcl::PointCloud<pcl::PointXYZ>::Ptr cyd(new pcl::PointCloud<pcl::PointXYZ>);
	//for ()
	////for (int j = 0; j < syd_candidate[i].size(); j++){ (*zzzzz).push_back((*cloud_src)[syd_candidate[i][j]]); }
	//visualize_pcd(zzzzz, cloud_no_plain_ptr, cloud_src);
    //筛选圆柱
	//首先从大到小排序
	for (int i = 0; i < syd_candidate.size()-1; i++){
		int maxids = i;
		for (int j = i + 1; j < syd_candidate.size(); j++)
		{
			maxids = syd_candidate[maxids] >= syd_candidate[j] ? maxids : j;
		}
		if (maxids==i) continue;
		int tmp = syd_candidate[maxids];
		syd_candidate[maxids] = syd_candidate[i];
		syd_candidate[i] = tmp;

		cylinder_coeffient_t * temp = cyd_cof_candi[maxids];
		cyd_cof_candi[maxids] = cyd_cof_candi[i];
		cyd_cof_candi[i] = temp;

		std::vector<float> ttemp = axis_length_candi[maxids];
		axis_length_candi[maxids] = axis_length_candi[i];
		axis_length_candi[i] = ttemp;
	}

	//合并圆柱
	//原因：实验过程中发现有些点数较少的干扰
	//合并圆柱 让正确的圆柱的点数更多

	std::vector<std::vector<int>> Same_syd_clu;//用来记录聚类的圆柱,记录的是在原数组中的索引
	std::vector<int>  Same_syd_clu_size ;//每一类圆柱总的点数
	std::vector<int> flag = std::vector<int>(syd_candidate.size(), 0);
	for (int i = 0; i < syd_candidate.size(); i++){
		if (flag[i]==1) continue;
		flag[i] = 1;
		std::vector<int> syd_clu_tmp;
		int size_tmp=0;
		syd_clu_tmp.push_back(i);
		size_tmp += syd_candidate[i];
		cylinder_coeffient_t * first_cyd = cyd_cof_candi[i];
		pcl::PointXYZ a1;

		a1.x = first_cyd->a;
		a1.y = first_cyd->b;
		a1.z = first_cyd->c;

		pcl::PointXYZ a2;
		a2.x = first_cyd->x;
		a2.y = first_cyd->y;
		a2.z = first_cyd->z;

		pcl::PointXYZ a3 = Pointadd(a2, Pointmulti(a1, axis_length_candi[0][0]));
		pcl::PointXYZ a4 = Pointadd(a2, Pointmulti(a1, axis_length_candi[0][1]));
		pcl::PointXYZ A, B, C;
		for (int j = i + 1; j < syd_candidate.size(); j++){
			if (flag[j]==1) continue;
			cylinder_coeffient_t * tmp = cyd_cof_candi[j];
			pcl::PointXYZ b1;
			b1.x = tmp->a;
			b1.y = tmp->b;
			b1.z = tmp->c;

			pcl::PointXYZ b2;
			b2.x = tmp->x;
			b2.y = tmp->y;
			b2.z = tmp->z;

			pcl::PointXYZ b3 = Pointadd(b2, Pointmulti(b1, axis_length_candi[j][0]));
			pcl::PointXYZ b4 = Pointadd(b2, Pointmulti(b1, axis_length_candi[j][1]));
			float angle_cos = PointDot(a1, b1) / (Pointlength(a1)*Pointlength(b1));
			std::vector <double> res;
			line2line_dis(a3, b3, a4, b4, res);

			//pcl::PointCloud<pcl::PointXYZ>::Ptr zzzzz(new pcl::PointCloud<pcl::PointXYZ>);
			//for (int j = 0; j < syd_candidate[i].size(); j++){ (*zzzzz).push_back((*cloud_src)[syd_candidate[i][j]]); }
			//visualize_pcd(cloud_no_plain_ptr, cloud_no_plain_ptr, zzzzz);
			if ((1-abs(angle_cos)) < 0.01&&abs(res[0]) < 0.1){
				flag[j] = 1;
				syd_clu_tmp.push_back(j);
				size_tmp += syd_candidate[j];
			}
		}
		Same_syd_clu.push_back(syd_clu_tmp);
		Same_syd_clu_size.push_back(size_tmp);
	}

    //按照每类的总点数排序
	std::vector<int> cyd_clu_index;
	for (int i = 0; i < Same_syd_clu_size.size(); i++) cyd_clu_index.push_back(i);
	for (int i = 0; i < Same_syd_clu_size.size()-1; i++){
		int min_size = i;
		for (int j = i + 1; j < Same_syd_clu_size.size(); j++)
		{
			min_size = Same_syd_clu_size[min_size] >= Same_syd_clu_size[j] ? min_size : j;
		}
		if (min_size == i) continue;

		int tmp_size = Same_syd_clu_size[min_size];
		Same_syd_clu_size[min_size] = Same_syd_clu_size[i];
		Same_syd_clu_size[i] = tmp_size;

		int tmp_index = cyd_clu_index[min_size];
		cyd_clu_index[min_size] = cyd_clu_index[i];
		cyd_clu_index[i] = tmp_index;
	}
    //挑出正确的圆柱
	pcl::PointCloud<pcl::PointXYZ>::Ptr center1(new pcl::PointCloud<pcl::PointXYZ>);
	//cyd_clu_index 类的size从大到小的索引 
	//Same_syd_clu  聚类之后的类 二维数组 ，每一行就是一个类，每行中存放的是syd_candidate的索引
	//由于最开始就是从到大小排序，所以每次取每一类的第一个就是该类点数最多的圆柱
	//默认最大的类的size最大的圆柱就一定是一个圆柱 然后通过是否垂直来找接下来的类的圆柱
	int index_first=Same_syd_clu[cyd_clu_index[0]][0];
	cylinder_coeffient_t * first_cyd = cyd_cof_candi[index_first];
	cylinder_coeffient_t * second_cyd=new cylinder_coeffient_t;
	pcl::PointXYZ centerPoint;
	pcl::PointXYZ a1;

	a1.x = first_cyd->a;
	a1.y = first_cyd->b;
	a1.z = first_cyd->c;
	
	pcl::PointXYZ a2;
	a2.x = first_cyd->x;
	a2.y = first_cyd->y;
	a2.z = first_cyd->z;

	pcl::PointXYZ a3 = Pointadd(a2, Pointmulti(a1, axis_length_candi[index_first][0]));
	pcl::PointXYZ a4 = Pointadd(a2, Pointmulti(a1, axis_length_candi[index_first][1]));
	pcl::PointXYZ A, B, C;
	int j=-1;
    for (int i = 1; i < cyd_clu_index.size(); i++){
		for (int k = 0; k < Same_syd_clu[cyd_clu_index[i]].size(); k++)
		{
			int index_second = Same_syd_clu[cyd_clu_index[i]][k];
			cylinder_coeffient_t * tmp = cyd_cof_candi[index_second];
			pcl::PointXYZ b1;
			b1.x = tmp->a;
			b1.y = tmp->b;
			b1.z = tmp->c;

			pcl::PointXYZ b2;
			b2.x = tmp->x;
			b2.y = tmp->y;
			b2.z = tmp->z;

			pcl::PointXYZ b3 = Pointadd(b2, Pointmulti(b1, axis_length_candi[index_second][0]));
			pcl::PointXYZ b4 = Pointadd(b2, Pointmulti(b1, axis_length_candi[index_second][1]));
			float angle_cos = PointDot(a1, b1) / (Pointlength(a1)*Pointlength(b1));
			std::vector <double> res;
			line2line_dis(a3, b3, a4, b4, res);

			if (abs(angle_cos) < 0.087&&abs(res[0]) < 0.2){
				(*center1).push_back(a3);
				(*center1).push_back(a4);
				(*center1).push_back(b3);
				(*center1).push_back(b4);
				second_cyd = tmp;
				centerPoint.x = res[1];
				centerPoint.y = res[2];
				centerPoint.z = res[3];
				j = index_second;
				i = cyd_clu_index.size();
				break;
			}
		}

	}
	if (j != -1)
    {
		cyd_vertical.push_back(first_cyd);
		cyd_vertical.push_back(second_cyd);
		pcl::Normal a;
		a.normal_x = first_cyd->a;
		a.normal_y = first_cyd->b;
		a.normal_z = first_cyd->c;
		pcl::Normal b;
		b.normal_x = second_cyd->a;
		b.normal_y = second_cyd->b;
		b.normal_z = second_cyd->c;
		pcl::Normal c;
		c.normal_x = -first_cyd->a;
		c.normal_y = -first_cyd->b;
		c.normal_z = -first_cyd->c;
		pcl::Normal d;
		d.normal_x = -second_cyd->a;
		d.normal_y = -second_cyd->b;
		d.normal_z = -second_cyd->c;
		Eigen::Affine3f transform1;
		Eigen::Affine3f transform2;
		Eigen::Affine3f transform3;
		Eigen::Affine3f transform4;
		twoaxis2trans(centerPoint, a, b, transform1);
		twoaxis2trans(centerPoint, b, c, transform2);
		twoaxis2trans(centerPoint, c, d, transform3);
		twoaxis2trans(centerPoint, d, a, transform4);
		pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud1(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud2(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud3(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_cloud4(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr tmp_cloud(new pcl::PointCloud<pcl::PointXYZ>);

		//(*center1).push_back(centerPoint);



		// 对点云进行变换
		pcl::transformPointCloud(*cloud_src, *transformed_cloud1, transform1);
		pcl::transformPointCloud(*cloud_src, *transformed_cloud2, transform2);
		pcl::transformPointCloud(*cloud_src, *transformed_cloud3, transform3);
		pcl::transformPointCloud(*cloud_src, *transformed_cloud4, transform4);

		pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_tgt(new pcl::search::KdTree<pcl::PointXYZ>());
		tree_tgt->setInputCloud(cloud_tgt);

		double err_dis1 = 0;
		double err_dis2 = 0;
		double err_dis3 = 0;
		double err_dis4 = 0;

		for (int i = 0; i < (*transformed_cloud1).size(); i++){
			std::vector<int> k_indices;
			std::vector<float> k_sqr_distances;
			tree_tgt->nearestKSearch((*transformed_cloud1)[i], 1, k_indices, k_sqr_distances);
			err_dis1 += k_sqr_distances[0];
		}

		for (int i = 0; i < (*transformed_cloud2).size(); i++){
			std::vector<int> k_indices;
			std::vector<float> k_sqr_distances;
			tree_tgt->nearestKSearch((*transformed_cloud2)[i], 1, k_indices, k_sqr_distances);
			err_dis2 += k_sqr_distances[0];
		}

		for (int i = 0; i < (*transformed_cloud3).size(); i++){
			std::vector<int> k_indices;
			std::vector<float> k_sqr_distances;
			tree_tgt->nearestKSearch((*transformed_cloud3)[i], 1, k_indices, k_sqr_distances);
			err_dis3 += k_sqr_distances[0];
		}

		for (int i = 0; i < (*transformed_cloud4).size(); i++){
			std::vector<int> k_indices;
			std::vector<float> k_sqr_distances;
			tree_tgt->nearestKSearch((*transformed_cloud4)[i], 1, k_indices, k_sqr_distances);
			err_dis4 += k_sqr_distances[0];
		}

		Eigen::Affine3f transform_res;
		double min_dis = std::min(std::min(err_dis1, err_dis2), std::min(err_dis3, err_dis4));
		if (min_dis == err_dis1) {
			for (int i = 0; i < transformed_cloud1->size();++i)
			cloud_final->push_back(transformed_cloud1->at(i));
			transform_res = transform1; }
		else if (min_dis == err_dis2) {
			for (int i = 0; i < transformed_cloud2->size(); ++i)
				cloud_final->push_back(transformed_cloud2->at(i));
			transform_res = transform2; }
		else if (min_dis == err_dis3){ 
			for (int i = 0; i < transformed_cloud3->size(); ++i)
				cloud_final->push_back(transformed_cloud3->at(i));
			transform_res = transform3; }
		else { 
			for (int i = 0; i < transformed_cloud4->size(); ++i)
				cloud_final->push_back(transformed_cloud4->at(i));
			transform_res = transform4; }

		Transform = std::vector<std::vector<float>>(4, std::vector<float>(4, 0));
		for (int row = 0; row < 4; row++){
			for (int col = 0; col < 4; col++){
				Transform[row][col] = transform_res(row, col);
			}
		}
		//Transform[0][3] = centerPoint.x;
		//Transform[1][3] = centerPoint.y;
		//Transform[2][3] = centerPoint.z;
		//Transform[3][3] = 1;

		////点云显示
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cyd_ptr(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_no_cyd_ptr(new pcl::PointCloud<pcl::PointXYZ>);
	}
    else{

		std::cout << "没有 垂直的 圆柱" << endl;
		//system("pause");
	}

}
