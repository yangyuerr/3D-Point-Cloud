
/*****************************************************************************
*                                 pcl_Eigen.c
*C语言实现ICP、下采样、法向量计算、pcd文件读取
*
*
*****************************************************************************/
#include <math.h>

//4*4矩阵相乘
void MatrixMul4D(double(*LMul)[4], double(*RMul)[4], double(*ResMat)[4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			ResMat[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				ResMat[i][j] += LMul[i][k] * RMul[k][j];
			}
		}
	}
}
/*
读取pcd文件
filename：文件路径
points：存储点云坐标到容器中
N：点云数量
*/
void readPCDfile(const char *filename, Vector *points, int *_N)
{
	FILE *fp;
	fopen_s(&fp,filename, "rb");
	if (fp == NULL) {
		printf("打开文件失败！");
		return;
	}

	char s[11][1024]; //存储Header
	int Points_Num; //点数量
	char data_columns_type[7] = { "" }; //列数: X Y Z
	char data_type[7] = { "" }; //点坐标存储方式(ascii或者binary)

	//连续读取Header行数   
	printf("start to read file header.....\n");
	printf("file header: \n");
	//以下11行是文件头信息
	for (int i = 0; i < 11; ++i) {
		fgets(s[i], 1024, fp);  //读取一行
		//printf("%s", s[i]); //输出

		//FIELDS x y z rgb
		if (i == 2) {
			char *s2 = s[2];
			char *pos = strstr(s2, "FIELDS");
			size_t size = strlen(s2);
			int n = min(size - 1 - 7, 7); // size - 1 - 7 > 7 ? size - 1 - 7 : 7;
			//strncpy(data_columns_type, pos + 7, n);
			//data_columns_type[n] = '\0';
			//printf("%s\n", data_columns_type);
		}

		//POINTS xxx
		if (i == 9) {
			char *s9 = s[9];
			char *pos = strstr(s9, "POINTS");
			size_t size = strlen(s9);
			int n = size - 1 - 7;
			char tmp[10] = { "" };
			strncpy(tmp, pos + 7, n);
			Points_Num = atoi(tmp);
			*_N = Points_Num;
			printf("Points_Num=%d\n", Points_Num);
		}

		//DATA ascii或者binary
		if (i == 10) {
			char *s10 = s[10];
			char *pos = strstr(s10, "DATA");
			size_t size = strlen(s10);
			int n = size - 1 - 5;
			strncpy(data_type, pos + 5, n);
			printf("data_type=%s\n", data_type);
		}

	}

	//以下开始读取xyz数据 
	//printf("start to read point .....\n");
	PointXYZ p;
	//vector_setup(&points, Points_Num, sizeof(PointXYZ));
	//if ((data_columns_type == "x y z") && (data_type == "ascii")) {
	//读取点坐标记录
	int ind = 0;
	if (1) {//(data_columns_type == "x y z") && (data_type == "ascii")
		while (!feof(fp)) {
			ind++;
			//if (ind>10000) {
			//	break;
			//}
			char tmp[1024];
			fgets(tmp, 1024, fp);  //读取一行
			p.x = atof(tmp);
			char *pos = strstr(tmp, " ");
			p.y = atof(pos);
			char *pos2 = strstr(pos + 1, " ");
			p.z = atof(pos2);
			vector_push_back(points, &p);
		}
	}
	else {
		//printf("data_type = binary, read failed!");
	}
}

/*
应用RT变换
points：原始点云
points_RT：变换后的点云
R：3*3旋转矩阵
T：1*3平移矩阵
*/
void applyRT(Vector points, Vector *points_RT, double R[][3], double T[], int _N) {
	PointXYZ p2;
	PointXYZ p_rt;
	for (size_t i = 0; i < _N; ++i)
	{
		p2 = *(PointXYZ*)vector_get(&points, i);
		float x, y, z;
		x = p2.x;
		y = p2.y;
		z = p2.z;

		p_rt.x = R[0][0] * x + R[0][1] * y + R[0][2] * z + T[0];
		p_rt.y = R[1][0] * x + R[1][1] * y + R[1][2] * z + T[1];
		p_rt.z = R[2][0] * x + R[2][1] * y + R[2][2] * z + T[2];
		vector_push_back(points_RT, &p_rt);
	}
}

/*
计算两点之间的距离
a1：原始点云中的点
a2：目标点云中的点
dim：点坐标维度
*/
static double dist_sq(double *a1, double *a2, int dims) {
	double dist_sq = 0, diff;
	while (--dims >= 0) {
		diff = (a1[dims] - a2[dims]);
		dist_sq += diff*diff;
	}
	return dist_sq;
}


/*----------------------------
*功能：找到输入点云中的包围盒两个点的值（右上和左下）
*-----------------------------
*输入：vector<PointXYZ>& InputCloudPoint（Piont3D的原始点云数据）
*输出：点云的min_p和max_p
*/
void GetMaxMin(Vector* InputCloudPoint, int size_input,PointXYZ*  min_p, PointXYZ* max_p)
{
	//主要思路是找到x，y,z的最小值,这样就能得到点云立体包围的尺寸
	//找x,y,z最小值
	if (size_input == 0)
	{
		printf("输入点云为空");
		return;
	}
	float x_min = INT_MAX, y_min= INT_MAX , z_min = INT_MAX;
	float x_max = INT_MIN, y_max = INT_MIN, z_max = INT_MIN;
	for (int i = 0; i < size_input; i++) 
	{
		PointXYZ a_tmp=*(PointXYZ*)vector_get(InputCloudPoint, i);
		x_min = a_tmp.x < x_min ? a_tmp.x : x_min;
		y_min = a_tmp.y < y_min ? a_tmp.y : y_min;
		z_min = a_tmp.z < z_min ? a_tmp.z : z_min;
		x_max = a_tmp.x > x_max ? a_tmp.x : x_max;
		y_max = a_tmp.y > y_max ? a_tmp.y : y_max;
		z_max = a_tmp.z > z_max ? a_tmp.z : z_max;
	}
	//给min_p赋值
	min_p->x = x_min;
	min_p->y = y_min;
	min_p->z = z_min;
	//给max_p赋值
	max_p->x = x_max;
	max_p->y = y_max;
	max_p->z = z_max;
	
	return;
}

/*----------------------------
*功能：体素化网格方法实现下采样（PCL中的源码C实现）
*-----------------------------
*输入：vector<PointXYZ>&InputCloudPoint (Piont3D的原始点云数据),下采样的体素大小X_Voxel,Y_Voxel,Z_Voxel
*输出：vector<PointXYZ>&OutPointCloud (采样之后的之后的PointXYZ结构的点云数据)
*/
void VoxelGrid_ApplyFilter(Vector* InputCloudPoint,int size_input,
	Vector* OutPointCloud, int* size_output,float X_Voxel, float Y_Voxel, float Z_Voxel)
{
	//先判断输入的点云是否为空
	if (size_input == 0)
	{
		printf("输入点云为空");
		return;
	}
	
	//存放输入点云的最大与最小坐标
	PointXYZ min_p, max_p;
	GetMaxMin(InputCloudPoint, size_input,&min_p, &max_p);

	PointXYZ inverse_leaf_size_;
	inverse_leaf_size_.x = 1 / X_Voxel;
	inverse_leaf_size_.y = 1 / Y_Voxel;
	inverse_leaf_size_.z = 1 / Z_Voxel;


	//计算最小和最大边界框值
	PointXYZ min_b_, max_b_, div_b_, divb_mul_;
	min_b_.x = (int)(floor(min_p.x * inverse_leaf_size_.x));
	max_b_.x = (int)(floor(max_p.x * inverse_leaf_size_.x));
	min_b_.y = (int)(floor(min_p.y * inverse_leaf_size_.y));
	max_b_.y = (int)(floor(max_p.y * inverse_leaf_size_.y));
	min_b_.z = (int)(floor(min_p.z * inverse_leaf_size_.z));
	max_b_.z = (int)(floor(max_p.z * inverse_leaf_size_.z));

	//计算沿所有轴所需的分割数
	div_b_.x = max_b_.x - min_b_.x + 1;
	div_b_.y = max_b_.y - min_b_.y + 1;
	div_b_.z = max_b_.z - min_b_.z + 1;

	//设置除法乘数
	divb_mul_.x = 1;
	divb_mul_.y = div_b_.x;
	divb_mul_.z = div_b_.x * div_b_.y;

	//用于计算idx和pointcloud索引的存储
	Vector index_vector;
	vector_setup(&index_vector, size_input, sizeof(Pair));

	//第一步：遍历所有点并将它们插入到具有计算idx的index_vector向量中;具有相同idx值的点将有助于产生CloudPoint的相同点
	for (unsigned int i = 0; i < size_input; i++)
	{
		PointXYZ p_tmp=*(PointXYZ*)vector_get(InputCloudPoint, i);
		int ijk0 = (int)(floor(p_tmp.x * inverse_leaf_size_.x) - (float)(min_b_.x));
		int ijk1 = (int)(floor(p_tmp.y * inverse_leaf_size_.y) - (float)(min_b_.y));
		int ijk2 = (int)(floor(p_tmp.z * inverse_leaf_size_.z) - (float)(min_b_.z));

		//计算质心叶索引
		int idx = ijk0 * divb_mul_.x + ijk1 * divb_mul_.y + ijk2 * divb_mul_.z;
		Pair a_tmp;
		a_tmp.first = idx;
		a_tmp.second = i;
		vector_push_back(&index_vector, &a_tmp);
	}
	//第二步：使用表示目标单元格的值作为索引对index_vector向量进行排序;实际上属于同一输出单元格的所有点都将彼此相邻
	//这里用选择排序
	for (int i = 0; i < size_input - 1; i++)
	{
		int minids = i;
		Pair a_tmp, b_tmp;
		for (int j = i + 1; j < size_input; j++)
		{	
			a_tmp  = *(Pair*)vector_get(&index_vector, minids);
			b_tmp  = *(Pair*)vector_get(&index_vector, j);
			minids = a_tmp.first <= b_tmp.first ? minids : j;
		}
		//swap
		a_tmp = *(Pair*)vector_get(&index_vector, minids);
		b_tmp = *(Pair*)vector_get(&index_vector, i);
		vector_assign(&index_vector, i, &a_tmp);
		vector_assign(&index_vector, minids, &b_tmp);
	}


	//第三步：计数输出单元格，我们需要跳过所有相同的，相邻的idx值
	unsigned int total = 0;
	unsigned int index = 0;
	unsigned int min_points_per_voxel_ = 0;
	//first_and_last_indices_vector [i]表示属于对应于第i个输出点的体素的index_vector中的第一个点的index_vector中的索引，以及不属于第一个点的索引
	Vector first_and_last_indices_vector;
	vector_setup(&first_and_last_indices_vector, size_input, sizeof(Pair));                             //分配内存空间

	while (index < size_input)
	{
		unsigned int i = index + 1;
		while (i < size_input && (*(Pair*)vector_get(&index_vector, i)).first == (*(Pair*)vector_get(&index_vector, index)).first)
			++i;
		if (i - index >= min_points_per_voxel_)
		{
			++total;
			Pair a_tmp;
			a_tmp.first = index;
			a_tmp.second = i;
			vector_push_back(&first_and_last_indices_vector, &a_tmp);
		}
		index = i;
	}

	//第四步：计算质心，将它们插入最终位置
	//OutPointCloud.resize(total);      //给输出点云分配内存空间
	float x_Sum, y_Sum, z_Sum;
	PointXYZ PointCloud;
	unsigned int first_index, last_index;
	for (unsigned int cp = 0; cp < total; ++cp)
	{
		// 计算质心 - 来自所有输入点的和值，这些值在index_vector数组中具有相同的idx值
		Pair a_tmp;
		a_tmp = *(Pair*)vector_get(&first_and_last_indices_vector, cp);
		first_index = a_tmp.first;
		last_index  = a_tmp.second;
		x_Sum = 0;
		y_Sum = 0;
		z_Sum = 0;
		for (unsigned int li = first_index; li < last_index; ++li)
		{
			Pair b_tmp;
			b_tmp = *(Pair*)vector_get(&index_vector, li);
			PointXYZ c_tmp;
			c_tmp = *(PointXYZ*)vector_get(InputCloudPoint, b_tmp.second);
			x_Sum += c_tmp.x;
			y_Sum += c_tmp.y;
			z_Sum += c_tmp.z;
		}
		PointCloud.x = x_Sum / (last_index - first_index);
		PointCloud.y = y_Sum / (last_index - first_index);
		PointCloud.z = z_Sum / (last_index - first_index);
		vector_push_back(OutPointCloud,&PointCloud);
	}
	*size_output = total;
	return;
}
double an_abs(double x) {
	if (x>0) {
		return x;
	}
	else {
		return -x;
	}
}

// ref: https://blog.csdn.net/album_gyd/article/details/81416398
double SqrtByNewton(double x)
{
	double val = x;//最终
	double eps = 0.00001;
	double last;//保存上一个计算的值
	do
	{
		last = val;
		val = (val + x / val) / 2;
	} while (an_abs(val - last) > eps);
	return val;
}

void vector_normalize(double* a, int len) {
	double sum = 0.0;
	for (int i = 0; i<len; i++) {
		sum += a[i] * a[i];
	}
	sum = SqrtByNewton(sum);
	for (int i = 0; i<len; i++) {
		a[i] /= sum;
	}
}

// if a == b, return 1
// if a != b, return 0
int vector_equal(double* a, double* b, int len) {
	int is_equal = 0;
	double eps = 0.001;

	for (int i = 0; i<len; i++) {
		int res = (an_abs(a[i] - b[i]) <= eps);
		is_equal += res;
	}

	return (is_equal == len);
}

// computation normal with two nearest points
int sub_normal_computation(double* tp, double* np1, double* np2, double* n) {
	int is_valid = -1;

	n[0] = 0;
	n[1] = 0;
	n[2] = 0;

	// it cannot estimate normal if tp=np1 or tp=np2 or np1=np2
	int s1 = vector_equal(tp, np1, 3);
	int s2 = vector_equal(tp, np2, 3);
	int s3 = vector_equal(np1, np2, 3);

	if (s1 + s2 + s3 > 0) {
		return is_valid;
	}

	double v1[3], v2[3];
	for (int i = 0; i<3; i++) {
		v1[i] = np1[i] - tp[i];
		v2[i] = np2[i] - tp[i];
	}

	vector_normalize(v1, 3);
	vector_normalize(v2, 3);
	// printf("\n there \n");
	// for(int a=0; a<3; a++){printf("%f ",v1[a]);}
	// for(int a=0; a<3; a++){printf("%f ",v2[a]);}
	// printf("\n --- \n");

	// it cannot estimate normal if vector tp-np1 conicide with vector tp-np2
	int s4 = vector_equal(v1, v2, 3);
	if (s4 > 0) {
		return is_valid;
	}

	is_valid = 1;

	// calculate normal now
	// n = v1*v2 where * is vector product
	n[0] = v1[1] * v2[2] - v1[2] * v2[1];
	n[1] = (v1[0] * v2[2] - v1[2] * v2[0])*(-1);
	n[2] = v1[0] * v2[1] - v1[1] * v2[0];
	vector_normalize(n, 3);

	if (n[2] <= 0) {
		n[0] *= (-1);
		n[1] *= (-1);
		n[2] *= (-1);
	}

	return is_valid;
}

// function for normal computation
// tp  is target point position, [1,3] vector
// np  is nearest points position, [N,3] matrix
// len is the number of nearest points
// n   is the normal result, [1,3] vector
int normal_computation(double* tp, double np[][3], const int len, double* n) {
	int is_valid = -1;

	//for (int k = 0; k < 3; k++){
	//	printf("tp %f  ,", tp[k]);
	//}
	//printf("\n");

	// step one: compute relative position of nearest points
	//double mat_A[len][3];
	for (int i = 0; i < len; i++){

		//for (int k = 0; k < 3; k++){
		//	printf("np %f  ,", np[i][k]);
		//}
		//printf("\n");

		np[i][0] = np[i][0] - tp[0];
		np[i][1] = np[i][1] - tp[1];
		np[i][2] = np[i][2] - tp[2];

		//for (int k = 0; k < 3; k++){
		//	printf("%f  ,", np[i][k]);
		//}
		//printf("\n");
	}

	double mat_ata[3][3]; // mat_ata = a.t() * a
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			mat_ata[i][j] = 0.0;
			for (int k = 0; k < len-1; k++){
				mat_ata[i][j] += np[k][i] * np[k][j];
				//printf("debug-bef: %f ", mat_ata[i][j]);
				//printf("\n");
			}
		}
	}

	// step two: get 3*3 matrix ATA
	double **mat_ATA;
	mat_ATA = fmatrix(1, 3, 1, 3); // note: index start is 1 not 0
	for (int i = 1; i <= 3; i++){
		for (int j = 1; j <= 3; j++){
			mat_ATA[i][j] = mat_ata[i - 1][j - 1];
	//		printf(" %f ", mat_ATA[i][j]);
		}
	//	printf("\n");
	}

	// step three: prepare for svd: MAT = u.t() * Diag(w) * v
	double wmax, wmin, **u, *w, **v;
	u = fmatrix(1, 3, 1, 3);
	v = fmatrix(1, 3, 1, 3);
	w = fvector(1, 3);

	//printf("u-bef : ");
	for (int i = 1; i <= 3; i++){
		for (int j = 1; j <= 3; j++){
			u[i][j] = mat_ATA[i][j];
		//	printf(" %f ", mat_ATA[i][j]);
		}
		//printf("\n");
	}

	// step four: svd
	svdcmp(u, 3, 3, w, v);

//	printf("u : ");
	for (int i = 1; i <= 3; i++){
		for (int j = 1; j <= 3; j++){
		//	printf(" %f ", u[i][j]);
		}
		//printf("\n");
	}

	//printf("v : ");
	for (int i = 1; i <= 3; i++){
		for (int j = 1; j <= 3; j++){
		//	printf(" %f ", v[i][j]);
		}
		//printf("\n");
	}


	//printf("w : ");
	//for (int i = 1; i <= 3; i++){ printf("%f   ", w[i]); }
	//printf("\n");


	int min_idx = -1;
	wmin = 1e7;
	for (int j = 1; j <= 3; j++){
		if (fabs(w[j]) < wmin){
			wmin = fabs(w[j]);
			min_idx = j;
		}
	}
	//printf("wmin: %d\n", wmin);
	//printf("w[0] %d\n", w[1]);
	//printf("w[1] %d\n", w[2]);
	//printf("w[2] %d\n", w[3]);
	//printf("%d\n", min_idx);
	// step five: find normal 
	for (int i = 1; i <= 3; i++){
		n[i - 1] = v[i][min_idx];
	}

	// step six: change direction if needed
	if (n[2] <= 0){
		n[0] *= (-1);
		n[1] *= (-1);
		n[2] *= (-1);
	}
	//printf("normasl : ");
	//for (int i = 0; i < 3; i++){ printf("%f   ", n[i]); }
	//printf("\n");
	//// reports
	//printf("eigenvalue: ");
	//for (int i = 1; i < 4; i++){ printf("%f ", w[i]); }
	//printf("\n");
	//printf("minimal eigenvalue: ");
	//printf("%f \n", wmin);
	//printf("minimal index: ");
	//printf("%d \n", min_idx);
	//printf("normal: ");
	//for (int i = 0; i < 3; i++){ printf("%f ", n[i]); }
	//printf("\n");

	return is_valid;
}

void Point_normal(Vector *PointCloud, int size, Vector *PointNormals, struct kdtree * kd,double radius)
{
	//PointXYZ pointA;
	//for (int i = 0; i < size; i++)
	//{
	//	pointA = *(PointXYZ*)vector_get(PointCloud, i);
	//	printf("x  %f", pointA.x);
	//	printf("y  %f", pointA.y);
	//	printf("z  %f", pointA.z);
	//	printf("\n");
	//}
	for (int i = 0; i < size; i++)
	{	
		struct kdres *presults;
		double pt[3],pos[3];
		double normal_target_point[3];
		PointNormal n_tmp;
		PointXYZ p_tmp = *(PointXYZ*)vector_get(PointCloud, i);
		pt[0] = p_tmp.x;
		pt[1] = p_tmp.y;
		pt[2] = p_tmp.z;
		presults=kd_nearest_range(kd, pt, radius);
		double nearest_points[1000][3];
		int j = 0;
		while (!kd_res_end(presults)) {
			
			//获取当前结果项的数据和位置
			(char*)kd_res_item(presults, pos);
			nearest_points[j][0] = pos[0];
			nearest_points[j][1] = pos[1];
			nearest_points[j][2] = pos[2];
			//printf("pos[0]   %f",pos[0]);
			//printf("pos[1]   %f", pos[1]);
			//printf("pos[2]   %f", pos[2]);
			//printf("\n");
			//转到下一个条目
			j++;
			//printf("j=%d\n", j);
			kd_res_next(presults);
		}
		normal_computation(pt, nearest_points, j, normal_target_point);
		n_tmp.normal_x = normal_target_point[0];
		n_tmp.normal_y = normal_target_point[1];
		n_tmp.normal_z = normal_target_point[2];
		vector_push_back(PointNormals, &n_tmp);

		kd_res_free(presults);
	}


}

/*icp算法
*输入原点云 原点云的索引，原点云索引的数量，目标点云及数量，迭代次数,以及变换方程
*输出平均点到点的距离误差
*/

double PointCloud_Icp(Vector *Cloud_src, Vector *index, int size_index, Vector *Cloud_tgt, int size_tgt, int iter, double (*Trans_icp)[4])
{
	//首先建立kdtree
	struct kdtree *kd_tgt;
	kd_tgt = kd_create(3);
	PointXYZ p_tmp;
	for (size_t i = 0; i < size_tgt; ++i)
	{
		p_tmp = *(PointXYZ*)vector_get(Cloud_tgt, i);
		assert(kd_insert3(kd_tgt, p_tmp.x, p_tmp.y, p_tmp.z, 0) == 0);
	}
	Vector Clould_src_noP;
	vector_setup(&Clould_src_noP, 10000, sizeof(PointXYZ));
	for (size_t i = 0; i < size_index; ++i)//points.size()
	{
		PointXYZ p_p;
		int index_tmp = *(int*)vector_get(index, i);
		p_p = *(PointXYZ*)vector_get(Cloud_src, index_tmp);

		vector_push_back(&Clould_src_noP, &p_p);
	}

	while (iter--)
	{
		//printf("iter=%d\n", iter);
		//然后建立最近点的方程
		double A[3][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		//计算中心点,为建立方程做准备
		double m_x = 0, m_y = 0, m_z = 0;
		double m_x_n = 0, m_y_n = 0, m_z_n = 0;
		PointXYZ p_p;
		PointXYZ p_n;
		Vector points_nearest;
		vector_setup(&points_nearest, 10000, sizeof(PointXYZ));
		for (size_t i = 0; i < size_index; ++i)//points.size()
		{
			p_p = *(PointXYZ*)vector_get(&Clould_src_noP, i);
			struct kdres *presult = kd_nearest(kd_tgt, &p_p);
			double pos[3];
			(char*)kd_res_item(presult, pos);
			m_x += p_p.x;
			m_y += p_p.y;
			m_z += p_p.z;
			m_x_n += pos[0];
			m_y_n += pos[1];
			m_z_n += pos[2];
			p_n.x = pos[0];
			p_n.y = pos[1];
			p_n.z = pos[2];
			vector_push_back(&points_nearest, &p_n);
			kd_res_free(presult);
		}
		m_x /= size_index;
		m_y /= size_index;
		m_z /= size_index;
		m_x_n /= size_index;
		m_y_n /= size_index;
		m_z_n /= size_index;



		for (size_t i = 0; i < size_index; ++i)//points.size()
		{
			p_p = *(PointXYZ*)vector_get(&Clould_src_noP, i);
			p_n = *(PointXYZ*)vector_get(&points_nearest, i);
			A[0][0] += (p_n.x - m_x_n)*(p_p.x - m_x);
			A[0][1] += (p_n.x - m_x_n)*(p_p.y - m_y);
			A[0][2] += (p_n.x - m_x_n)*(p_p.z - m_z);
			A[1][0] += (p_n.y - m_y_n)*(p_p.x - m_x);
			A[1][1] += (p_n.y - m_y_n)*(p_p.y - m_y);
			A[1][2] += (p_n.y - m_y_n)*(p_p.z - m_z);
			A[2][0] += (p_n.z - m_z_n)*(p_p.x - m_x);
			A[2][1] += (p_n.z - m_z_n)*(p_p.y - m_y);
			A[2][2] += (p_n.z - m_z_n)*(p_p.z - m_z);

		}
		//p_p = *(PointXYZ*)vector_get(&Clould_src_noP, 0);
		//p_n = *(PointXYZ*)vector_get(&points_nearest, 0);
		//printf("%f ", p_p.x);
		//printf("%f ", p_p.y);
		//printf("%f ", p_p.z);
		//printf("%f ", p_n.x);
		//printf("%f ", p_n.y);
		//printf("%f ", p_n.z);
		//printf("\n");

		double U[3][3];
		double S[3];
		double V[3][3];
		int k = 3;
		svd(A, k, U, S, V);

		//UV
		double U_tmp[3][3], V_tmp[3][3];
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				U_tmp[i][j] = U[i][j];
				V_tmp[i][j] = V[i][j];
			}
		}

		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				U[i][j] = U_tmp[j][i];
				V[i][j] = V_tmp[j][i];
			}
		}

		U[0][0] = -U[0][0];
		U[0][2] = -U[0][2];
		U[1][0] = -U[1][0];
		U[1][2] = -U[1][2];
		U[2][0] = -U[2][0];
		U[2][2] = -U[2][2];

		V[0][0] = -V[0][0];
		V[0][2] = -V[0][2];
		V[1][0] = -V[1][0];
		V[1][2] = -V[1][2];
		V[2][0] = -V[2][0];
		V[2][2] = -V[2][2];


		//-----------------------计算R、T-----------------------//
		double R[3][3] = { 0,0,0,0,0,0,0,0,0 };
		double T[3] = { 0,0,0 };

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				R[i][j] = U[i][0] * V[j][0] + U[i][1] * V[j][1] + U[i][2] * V[j][2];
			}
		}

		T[0] = m_x_n - (R[0][0] * m_x + R[0][1] * m_y + R[0][2] * m_z);
		T[1] = m_y_n - (R[1][0] * m_x + R[1][1] * m_y + R[1][2] * m_z);
		T[2] = m_z_n - (R[2][0] * m_x + R[2][1] * m_y + R[2][2] * m_z);




		//for (int i = 0; i < 2; i++)
		//{
		//	for (int j = i + 1; j < 3; j++)
		//	{
		//		double tmp = R[i][j];
		//		R[i][j] = R[j][i];
		//		R[j][i] = tmp;
		//	}
		//}
		//double T_tmp[3];
		//T_tmp[0] = R[0][0] * T[0] + R[0][1] * T[1] + R[0][2] * T[2];
		//T_tmp[1] = R[1][0] * T[0] + R[1][1] * T[1] + R[1][2] * T[2];
		//T_tmp[2] = R[2][0] * T[0] + R[2][1] * T[1] + R[2][2] * T[2];

		//T[0] = T_tmp[0];
		//T[1] = T_tmp[1];
		//T[2] = T_tmp[2];
		//点云变换
		Vector points_tmp;
		vector_setup(&points_tmp, 10000, sizeof(PointXYZ));
		applyRT(Clould_src_noP, &points_tmp, R, T, size_index);
		vector_move_assign(&Clould_src_noP, &points_tmp);

		double Trans[4][4] = {R[0][0],R[0][1], R[0][2], T[0],
							  R[1][0], R[1][1], R[1][2],T[1] ,
							  R[2][0], R[2][1], R[2][2],T[2] ,
							  0,0,0,1};
		//变换矩阵
		double res[4][4];
		MatrixMul4D( Trans, Trans_icp,res);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Trans_icp[i][j] = res[i][j];
			}
		}
		
		//内存保护
		vector_clear(&points_nearest);
		vector_destroy(&points_nearest);
	}
	//计算均差
	PointXYZ p_p;
	double m_error=0;
	for (size_t i = 0; i < size_index; ++i)//points.size()
	{
		p_p = *(PointXYZ*)vector_get(&Clould_src_noP, i);
		struct kdres *presult = kd_nearest(kd_tgt, &p_p);
		double pos[3];
		(char*)kd_res_item(presult, pos);
		m_error += sqrt(pow((p_p.x - pos[0]), 2) + pow((p_p.y - pos[1]), 2) + pow((p_p.z - pos[2]), 2));
		kd_res_free(presult);
	}
	m_error /= size_index;

	//内存保护
	vector_clear(&Clould_src_noP);
	vector_destroy(&Clould_src_noP);
	kd_free(kd_tgt);

	return m_error;
}