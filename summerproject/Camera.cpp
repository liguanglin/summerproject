#include "Camera.h"
#include<iostream>

void Camera::Read(FILE* file)
{
	R = cv::Mat::zeros(3, 3, CV_32FC1);
	T = cv::Mat::zeros(3, 1, CV_32FC1);
	/*for (int i = 0; i < 3; i++) {
		double t; fscanf(file, "%lf", &t);
		T.at<float>(i, 0) = t;
	}
	double x, y, z, w;
	fscanf(file, "%lf%lf%lf%lf", &x, &y, &z, &w);
	Eigen::Quaternionf q(w, x, y, z);
	q.normalize();
	std::cout << q.w()<<" "<<q.x()<<" "<<q.y()<<" "<<q.z() << std::endl;
	cv::eigen2cv(q.matrix(), R);
	std::cout << R << std::endl;*/
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			double t; fscanf(file, "%lf", &t);
			if (j < 3) {
				R.at<float>(i, j) = t;
			}
			else {
				T.at<float>(i, 0) = t;
			}
		}
	}
	//std::cout << R << std::endl << T << std::endl;
}