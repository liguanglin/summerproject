#pragma once
#include<core/core.hpp>
#include<highgui/highgui.hpp> 
#include "opencv2/imgproc.hpp"

#include "Camera.h"
#include<iostream>

using namespace std;

class Frame {
public:
	cv::Mat color_img;
	cv::Mat gray_img;
	cv::Mat dep_img;
	int out;
	Camera *cam;
	int width, height;
	Frame() {}
	void setimg(cv::Mat&img);
	void readimg(char*file);
	void setcam(Camera*_cam);
	void GetWorldCoordFrmImgCoord(int x, int y, double d,cv::Mat &w_c);
	void GetImgCoordFrmWorldCoord(double&x, double&y, double&d, cv::Mat w_c);
};