#include "Camera.h"
#include<iostream>

void Camera::Read(FILE* file)
{
	R = cv::Mat::zeros(3, 3, CV_32FC1);
	T = cv::Mat::zeros(3, 1, CV_32FC1);
	for (int i = 0; i < 3; i++) {
		double t; fscanf(file, "%lf", &t);
		T.at<float>(i, 0) = t;
	}
	double x, y, z, w;
	fscanf(file, "%lf%lf%lf%lf", &x, &y, &z, &w);
	Eigen::Quaternionf q(w, x, y, z);
	q.normalize();
	std::cout << q.w()<<" "<<q.x()<<" "<<q.y()<<" "<<q.z() << std::endl;
	cv::eigen2cv(q.matrix(), R);
	std::cout << R << std::endl;
}