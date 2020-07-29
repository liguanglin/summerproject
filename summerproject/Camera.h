#pragma once
#include<core/core.hpp>
#include<highgui/highgui.hpp>

class Camera{
public:
	Camera() {};
	Camera(float _fx, float _fy, float _cx, float _cy, float _s=1.0) {
		in = cv::Mat::zeros(3, 3, CV_32FC1);
		in.at<float>(0, 0) = _fx;
		in.at<float>(1, 1) = _fy;
		in.at<float>(0, 2) = _cx;
		in.at<float>(1, 2) = _cy;
		in.at<float>(2, 2) = _s;
	}
	void Read(FILE*file);
	cv::Mat R,T,in;
};