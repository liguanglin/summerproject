#include<iostream>
#include<core/core.hpp>
#include<highgui/highgui.hpp>  
#include"Camera.h"
#include"Frame.h"
#include <Eigen\Dense>
#include "GCoptimization.h"

using namespace std;


int img_width, img_height;
double delta = 0.0001;
int k = 10;

int src(int x, int y, int label)
{
	return (x * img_width + y)*k + label;
}

void bp(double * &L_init,int* &result)
{
	/*double *nL = new double[img_width*img_height*k];
	double *x = L_init, *y = nL;
	int iter = 5;
	while (iter--) {
		for (int i = 0; i < img_height; i++) {
			for (int j = 0; j < img_width; j++) {
				for (int dy = 0; dy < k; dy++) {
					y[src(i, j, dy)] = x[src(i, j, dy)];
					for (int dx = 0; dx <= k; dx++) {
						y[src(i, j, dy)] += 1.0 / k * (x[src(i, j, dx)] + abs(dx*delta - dy * delta));
					}
				}
			}
		}
		swap(x, y);
	}
	L_init = x;
	delete y;*/
	int num_pixels = img_height * img_width;
	result = new int[num_pixels];   // stores result of optimization
	GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(img_width, img_height, k);
	for (int i = 0; i < img_height*img_width; i++) {
		for (int j = 0; j < k; j++) {
			gc->setDataCost(i, j, L_init[i*k+j]);
		}
	}
	for (int l1 = 0; l1 < k; l1++)
		for (int l2 = 0; l2 < k; l2++) {
			int cost = (l1 - l2)*(l1 - l2);
			gc->setSmoothCost(l1, l2, 0);
		}
	cout << gc->compute_energy();
	gc->expansion(5);
	cout << gc->compute_energy();
	for (int i = 0; i < num_pixels; i++)
		result[i] = gc->whatLabel(i);
}

void init(Frame*frames)
{
	int height =frames[0].height , width =frames[0].width;
	cout << width << " " << height << endl;
	img_width = width;
	img_height = height;
	double*L_init = new double[width*height*k];
	for (int i = 0; i < width*height*k; i++) {
		L_init[i] = 0;
	}
	cv::Mat trans=cv::Mat(height, width, CV_8UC1);
	for (int i = 0; i <= 0; i++) {
		cout << i << endl;
		Frame *frame = &frames[i];
		double ux = 1;
		for (int d = 1; d <= k; d++) {
			for (int j = 1; j < 5; j++) {
				Frame *framet = &frames[j];
				cv::Mat w_c0, w_c1, w_c2, w_c3;
				float x0, y0, x1, y1, x2, y2, x3, y3, tmpd;
				frame->GetWorldCoordFrmImgCoord(0, 0, 1 / (d*delta), w_c0);
				frame->GetWorldCoordFrmImgCoord(height, 0, 1 / (d*delta), w_c1);
				frame->GetWorldCoordFrmImgCoord(0, width, 1 / (d*delta), w_c2);
				frame->GetWorldCoordFrmImgCoord(height, width, 1 / (d*delta), w_c3);
				framet->GetImgCoordFrmWorldCoord(x0, y0, tmpd, w_c0);
				framet->GetImgCoordFrmWorldCoord(x1, y1, tmpd, w_c1);
				framet->GetImgCoordFrmWorldCoord(x2, y2, tmpd, w_c2);
				framet->GetImgCoordFrmWorldCoord(x3, y3, tmpd, w_c3);
				cout << x0 << " " << y0 << " \n" << x1 << " " << y1 << "\n" << x2 << " " << y2 <<"\n"<< x3 << " " << y3 << endl;
				//cout << tmpd << endl;
				for (int x = 0; x < height; x++) {
					for (int y = 0; y < width; y++) {
						double tmp = 0;
						int tx = (int)(x0 + 1.0*x / height * (x1 - x0)+0.5), ty = (int)(y0 + 1.0*y / width * (y2 - y0)+0.5);
						//std::cout << x << " " << y << " " << tx << " " << ty << endl;
						if (tx >= 0 && tx < height&&ty >= 0 && ty < width) {
							trans.at<uchar>(tx, ty) = frame->gray_img.at<uchar>(x, y);
							tmp = abs(framet->gray_img.at<uchar>(tx, ty) - frame->gray_img.at<uchar>(x, y));
							tmp = 1 / (1 + tmp * tmp);
						}

						L_init[src(x, y, d - 1)] += tmp;
						if (L_init[src(x, y, d - 1)] > ux)
							ux = L_init[src(x, y, d - 1)];
					}
				}
			}
			char file_name[20];
			sprintf(file_name, "%d.bmp", d);
			imwrite(file_name, trans);
		}
		ux = 100.0 / ux;
		//printf("%lf\n", ux);
		for (int x = 0; x < width*height*k; x++) {
			L_init[x] = 100 - ux * L_init[x];
			//cout << L_init[x]<<" ";
		}
		//cout << "bp begin\n" << endl;
		int *res;
		bp(L_init,res);
		//cout << "bp done\n" << endl;
		for (int x = 0; x < height; x++) {
			//cout << x << endl;
			for (int y = 0; y < width; y++) {
				//if(mind!=0) cout << x << " " << y << " " << mind * delta << endl;
				frame->dep_img.at<uchar>(x, y)= (res[x*width+y])*255/k;
			}
		}
		cv::imwrite("dep.bmp",frame->dep_img);
	}
	delete L_init;
}

int main()
{
	FILE* actfile = fopen("DATA/myact.txt", "r");
	int picnums;
	fscanf(actfile, "%d", &picnums);
	Frame* frames = new Frame[picnums];
	double fx,fy,cx,cy,s;
	fscanf(actfile,"%lf%lf%lf%lf%lf", &fx, &fy, &cx, &cy, &s);
	for (int i = 0; i < picnums; i++) {
		char filename[100];
		//fscanf(actfile,"%lf", &s);
		fscanf(actfile,"%s", filename);
		Camera*cam = new Camera(fx, fy, cx, cy);
		fscanf(actfile,"%lf", &s);
		cam->Read(actfile);
		frames[i].cam = cam;
		sprintf(filename, "DATA/test%04d.jpg", i);
		frames[i].readimg(filename);
		fscanf(actfile, "%lf,%lf,%lf,%lf", &s, &s, &s, &s);
		fscanf(actfile, "%s", filename);
	}
	init(frames);
}