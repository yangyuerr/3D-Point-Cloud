/*****************************************************************************
*                                rough_fine.cpp
*该文件主要是简单的点云粗、精匹配:
*为了方便，本文件只写核心代码，不封装
*
*此文件没有主函数，可在registration.cpp中调用此文件中的函数
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/

//NDT配准
//初始化正态分布变换（NDT）
pcl::NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> ndt;
//为终止条件设置最小转换差异
ndt.setTransformationEpsilon(0.05);
//为More-Thuente线搜索设置最大步长
ndt.setStepSize(0.07);
//设置NDT网格结构的分辨率（VoxelGridCovariance）（体素格的大小）
ndt.setResolution(0.7);
//设置匹配迭代的最大次数
ndt.setMaximumIterations(40);
// 设置要配准的点云
ndt.setInputSource(cloud_src);
//设置点云配准目标
ndt.setInputTarget(cloud_tgt_o);
//计算需要的刚体变换以便将输入的点云匹配到目标点云
pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
ndt.align(*output_cloud, sac_trans);



//SAC配准
pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
scia.setInputSource(cloud_src);
scia.setInputTarget(cloud_tgt);
scia.setSourceFeatures(fpfhs_src);
scia.setTargetFeatures(fpfhs_tgt);
//scia.setMinSampleDistance(1);
//scia.setNumberOfSamples(2);
//scia.setCorrespondenceRandomness(20);
PointCloud::Ptr sac_result(new PointCloud);
scia.align(*sac_result);
std::cout << "sac has converged:" << scia.hasConverged() << "  score: " << scia.getFitnessScore() << endl;
Eigen::Matrix4f sac_trans;
sac_trans = scia.getFinalTransformation();
std::cout << sac_trans << endl;
pcl::io::savePCDFileASCII("junjie_transformed_sac.pcd", *sac_result);



//icp配准
PointCloud::Ptr icp_result(new PointCloud);
pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
icp.setInputSource(cloud_src);
icp.setInputTarget(cloud_tgt_o);
//Set the max correspondence distance to 4cm (e.g., correspondences with higher distances will be ignored)
icp.setMaxCorrespondenceDistance(0.04);
// 最大迭代次数
icp.setMaximumIterations(50);
// 两次变化矩阵之间的差值
icp.setTransformationEpsilon(1e-10);
// 均方误差
icp.setEuclideanFitnessEpsilon(0.2);
icp.align(*icp_result, sac_trans);



//GICP
//PointXYZRGBNormal
// 成员变量:float x, y, z, rgb, normal[3], curvature;
// PointXYZRGBNormal存储XYZ数据和RGB颜色的point结构体，并且包括曲面法线和曲率。

// create the object implementing ICP algorithm
pcl::GeneralizedIterativeClosestPoint<pcl::PointXYZRGBNormal, pcl::PointXYZRGBNormal> gicp;
// set the input point cloud to align
gicp.setInputCloud(cloud_in);
// set the input reference point cloud
gicp.setInputTarget(cloud_out);
// compte the point cloud registration
pcl::PointCloud<pcl::PointXYZRGBNormal> Final;
gicp.align(Final);
// print if it the algorithm converged and its fitness score
std::cout << "has converged:" << gicp.hasConverged()
 << " score: "
 << gicp.getFitnessScore() << std::endl;
// print the output transformation
std::cout << gicp.getFinalTransformation() << std::endl;



