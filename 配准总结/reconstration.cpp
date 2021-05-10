/*****************************************************************************
*                                reconstruction.cpp
*该文件主要是简单的点云三角化重建：
*输入：pcd文件
*输出：三角化之后的ply文件
*作者：yangyuer   20210414  yangyu.hust@foxmail.com
*
*****************************************************************************/
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <boost/thread/thread.hpp>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
 
typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
int main(int argc, char** argv)
{
	//打开文件
	PointCloudT::Ptr cloud(new PointCloudT);
	if (pcl::io::loadPCDFile<pcl::PointXYZ>("face2.pcd", *cloud) == -1)
	{
		PCL_ERROR("Couldn't read file1 \n");
		return (-1);
	}
	// 估计法向量
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(cloud);
	n.setInputCloud(cloud);
	n.setSearchMethod(tree);
	n.setKSearch(20);
	n.compute(*normals); //计算法线，结果存储在normals中
	//* normals 不能同时包含点的法向量和表面的曲率
 
	//将点云和法线放到一起
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
	pcl::concatenateFields(*cloud, *normals, *cloud_with_normals);
	//* cloud_with_normals = cloud + normals
 
	//创建搜索树
	pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
	tree2->setInputCloud(cloud_with_normals);
 
	//初始化GreedyProjectionTriangulation对象，并设置参数
	pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
	//创建多变形网格，用于存储结果
	pcl::PolygonMesh triangles;
 
	//设置GreedyProjectionTriangulation对象的参数
	//第一个参数影响很大
	gp3.setSearchRadius(15.0f); //设置连接点之间的最大距离（最大边长）用于确定k近邻的球半径【默认值 0】
	gp3.setMu(5.0f); //设置最近邻距离的乘子，以得到每个点的最终搜索半径【默认值 0】
	gp3.setMaximumNearestNeighbors(100); //设置搜索的最近邻点的最大数量
	gp3.setMaximumSurfaceAngle(M_PI / 4); // 45 degrees（pi）最大平面角
	gp3.setMinimumAngle(M_PI / 18); // 10 degrees 每个三角的最小角度
	gp3.setMaximumAngle(2 * M_PI / 3); // 120 degrees 每个三角的最大角度
	gp3.setNormalConsistency(false); //如果法向量一致，设置为true
 
	//设置搜索方法和输入点云
	gp3.setInputCloud(cloud_with_normals);
	gp3.setSearchMethod(tree2);
 
	//执行重构，结果保存在triangles中
	gp3.reconstruct(triangles);
 
	//保存网格图
	pcl::io::savePLYFile("result.ply", triangles);
 
	// Additional vertex information
	//std::vector<int> parts = gp3.getPartIDs();
	//std::vector<int> states = gp3.getPointStates();
 
	// 显示结果图
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0); //设置背景
	viewer->addPolygonMesh(triangles, "my"); //设置显示的网格
	viewer->addCoordinateSystem(1.0); //设置坐标系
	viewer->initCameraParameters();
	while (!viewer->wasStopped()){
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
 
	return (0);
}