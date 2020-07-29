#include "Frame.h"

void Frame::setimg(cv::Mat&img)
{
	color_img = img.clone();
	width = img.cols;
	height = img.rows;
	gray_img = img.clone();
	cvtColor(color_img, gray_img, CV_BGR2GRAY);
	dep_img = cv::Mat(img.size(), CV_8UC1);
}
void Frame::readimg(char*file)
{
	color_img = cv::imread(file,CV_8UC1);
	width = color_img.cols;
	height = color_img.rows;
	gray_img = color_img.clone();
	//cvtColor(color_img,gray_img, CV_BGR2GRAY);
	dep_img = cv::Mat(color_img.size(), CV_8UC1);
}
void Frame::setcam(Camera*_cam)
{
	cam = _cam;
}

void Frame::GetWorldCoordFrmImgCoord(int x, int y, double d, cv::Mat &w_c)
{
	w_c = cv::Mat::zeros(3, 1, CV_32FC1);
	w_c.at<float>(0, 0) = x;
	w_c.at<float>(1, 0) = y;
	w_c.at<float>(2, 0) = 1;
	//cout << "------------------------" << endl;
	//cout << w_c << endl;
	w_c = d * cam->in.inv()*w_c;
	//cout << w_c << endl;
	//cout << cam->in.inv() << endl;
	w_c = cam->R.t()*w_c - cam->R.t()*cam->T;
	//cout << w_c << endl;
}

void Frame::GetImgCoordFrmWorldCoord(double&x, double&y, double&d, cv::Mat w_c)
{
	cv::Mat res;
	//cout << "++++++++++++++++++++++++" << endl;
	//cout << w_c << endl;
	res = cam->R*w_c + cam->T;
	//if(out)cout << res << endl;
	res = cam->in*res;
	//if(out)cout << res << endl;
	d = res.at<float>(2, 0);
	x = res.at<float>(0, 0)/d;
	y = res.at<float>(1, 0)/d;
}
