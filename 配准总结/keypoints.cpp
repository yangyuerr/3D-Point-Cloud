/*****************************************************************************
*                                keypoints.cpp
*该文件主要是简单的点云关键点提取:均匀采样、ISS、Harris3D、SIFT
*此文件没有主函数，可在registration.cpp中调用此文件中的函数
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/

//均匀采样关键点
//#include<pcl/filters/uniform_sampling.h>
computeUniformKeypoints(pcl::PointCloud<PointT>::Ptr input_src_cloud, pcl::PointCloud<PointT>::Ptr uniform_keypoint, double radius)
{
	pcl::UniformSampling<pcl::PointXYZ> uniform_sampler;
	uniform_sampler.setInputCloud(input_src_cloud);
	uniform_sampler.setRadiusSearch(radius);
	uniform_sampler.filter(*uniform_keypoint);
}

//ISS关键点
//#include<pcl/keypoints/iss_3d.h>
void CcomputeISSKeypoint(pcl::PointCloud<PointT>::Ptr input_cloud_ptr, pcl::PointCloud<PointT>::Ptr output_ISS_keypoint, double model_solution, int MinNeighbors)
{
	//double model_solution_tgt = 0.04;参数小，采取的关键点多，论文中为500左右
	//pcl::PointCloud<pcl::PointXYZ>::Ptr model_keypoint(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::ISSKeypoint3D<pcl::PointXYZ, pcl::PointXYZ> iss;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr iss_tree(new pcl::search::KdTree<pcl::PointXYZ>());
	//参数设置
	iss.setSearchMethod(iss_tree);
	iss.setSalientRadius(3 * model_solution);
	iss.setNonMaxRadius(2 * model_solution);
	iss.setThreshold21(0.975);
	iss.setThreshold32(0.975);
	iss.setMinNeighbors(MinNeighbors);
	iss.setNumberOfThreads(4);
	iss.setInputCloud(input_cloud_ptr);
	iss.compute(*output_ISS_keypoint);
}

//harris3D关键点
//#include<pcl/keypoints/harris_3d.h>
void computeHarrisKeypoints(pcl::PointCloud<PointT>::Ptr input_src_cloud, pcl::PointCloud<PointT>::Ptr harris_keypoint, double radius)
{
	pcl::HarrisKeypoint3D<pcl::PointXYZ, pcl::PointXYZI, pcl::Normal> harris_detector;
	harris_detector.setInputCloud(input_src_cloud);//设置输入点云 指针
	harris_detector.setNonMaxSupression(true);
	harris_detector.setRadius(radius);//　块体半径
	harris_detector.setThreshold(0.00000001f);//数量阈值
	//新建的点云必须初始化，清零，否则指针会越界
	//注意Harris的输出点云必须是有强度(I)信息的 pcl::PointXYZI，因为评估值保存在I分量里
	pcl::PointCloud<pcl::PointXYZI>::Ptr harris_keypointI(new pcl::PointCloud<pcl::PointXYZI>);

	// 计算特征点
	harris_detector.compute(*harris_keypointI);
	
	// 关键点
	pcl::copyPointCloud(*harris_keypointI, *harris_keypoint);
}

//sift关键点
//#include<pcl/keypoints/sift_keypoints.h>
namespace pcl
{
    template<>
    struct SIFTKeypointFieldSelector<PointXYZ>
    {
        /* data */
        inline float
            operator() (const PointXYZ &p) const
            {
                return p.z;
            }
    };
    
}
Void computeSIFT3DKeypoints()(pcl::PointCloud<PointT>::Ptr input_src_cloud,float min_scale, int n_octaves, int n_scales_per_octave, float min_contrast,pcl::PointCloud<PointT>::Ptr SIFT_keypoint)
{
    // Parameters for sift computation
    //const float min_scale = 0.0001f; //the standard deviation of the smallest scale in the scale space //越小点越多
    // const int n_octaves = 20;/,the number of octaves (i.e. doublings of scale) to compute
    //  const int n_scales_per_octave = 4;//the number of scales to compute within each octave //越小点越多
    //const float min_contrast = 0.0001f;//the minimum contrast required for detection

    //  Estimate the sift interest points using z values from xvz as the Intensity variants
    pcl::SIFTKeypoint<pcl::PointXYZ, pcl::PointWithScale> sift;
    pcl: :PointCloud<pcl: :PointWithScale> result;
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ> ());
    sift. setSearchMethod(tree);
    sift. setScales(min_scale, n_octaves, n_scales_per_octave);
    sift. setMinimumContrast(min_contrast);
    sift. setInputCloud(input_src_cloud);
    sift, compute(result);
    
    // Copying the ?oinmthscaLe to gointxvz so as visualize the cloud
    copyPointCloud(result, SIFT_keypoint);
}