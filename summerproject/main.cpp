#include<iostream>
#include<core/core.hpp>
#include<highgui/highgui.hpp>  
#include"Camera.h"
#include"Frame.h"
#include "GCoptimization.h"

using namespace std;


int img_width, img_height;
double delta = 0.001;
int k = 10;

int src(int x, int y, int label)
{
	return (x * img_width + y)*k + label;
}

struct ForDataFn {
	int numLab;
	int *data;
};

int dataFn(int p, int l, void *data)
{
	
	ForDataFn *myData = (ForDataFn *)data;
	int numLab = myData->numLab;
	return(myData->data[p*numLab + l]);
}

int smoothFn(int p1, int p2, int l1, int l2)
{
	if ((l1 - l2)*(l1 - l2) <= 100) return((l1 - l2)*(l1 - l2));
	else return(100);
}

void gc(int * data,int* &result)
{
	try {
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(img_width, img_height, k);
		ForDataFn toFn;
		printf("%d\n", data[5174400]);
		toFn.data = data;
		toFn.numLab = k;
		gc->setDataCost(&dataFn,&toFn);
		gc->setSmoothCost(&smoothFn);
		printf("\nBefore optimization energy is %lld", gc->compute_energy());
		gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
		printf("\nAfter optimization energy is %lld", gc->compute_energy());

		for (int i = 0; i < img_width*img_height; i++)
			result[i] = gc->whatLabel(i);
		delete gc;
	}
	catch (GCException e) {
		e.Report();
		while (1);
	}
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
	delete y;
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
		result[i] = gc->whatLabel(i);*/
}

double L_init[1000][1000][20];
double ux[1000][1000],u2[1000][1000];
double be[1000][1000][20];
double mess[2][1000][1000][20][4];
double summess[2][1000][1000][20];

double sqr(double x)
{
	return x * x;
}

double calsum(int x, int y, int x2, int y2, int iter, int d1, int d2,cv::Mat &frame,int dir)
{
	int It = frame.at<uchar>(x, y) - frame.at<uchar>(x2, y2);
	double sum = L_init[x2][y2][d2] + abs(d1 - d2)*u2[x][y] / abs(It+1e-7);
	sum += summess[iter ^ 1][x2][y2][d2] - mess[iter ^ 1][x2][y2][d2][dir];
	return sum;
}

void bp(cv::Mat &frame)
{
	memset(be, 0, sizeof(be));
	memset(summess, 0, sizeof(summess));
	memset(mess, 0, sizeof(mess));
	int n = img_height, m = img_width;
	int iter = 2;
	for (int iter = 0; iter < 2; iter++) {
		int ths=iter&1,las=(iter&1)^1;
		printf("%d\n", iter);
		for (int i = 1; i < n-1; i++) {
			//printf("%d\n", i);
			for (int j = 1; j < m-1; j++) {
				for (int d1 = 0; d1 < k; d1++) {
					double sum1 = 1e5, sum2 = 1e5, sum3 = 1e5, sum4 = 1e5;
					for (int d2 = 0; d2 < k; d2++) {
						sum1 = min(calsum(i, j, i + 1, j, ths, d1, d2, frame, 0), sum1);
						sum2 = min(calsum(i, j, i - 1, j, ths, d1, d2, frame, 1), sum2);
						sum3 = min(calsum(i, j, i, j - 1, ths, d1, d2, frame, 2), sum3);
						sum4 = min(calsum(i, j, i, j + 1, ths, d1, d2, frame, 3), sum4);
					}
					mess[ths][i][j][d1][0] = sum1;
					mess[ths][i][j][d1][1] = sum2;
					mess[ths][i][j][d1][2] = sum3;
					mess[ths][i][j][d1][3] = sum4;
					summess[ths][i][j][d1] = sum1 + sum2 + sum3 + sum4;
					//if(summess[ths][i][j][d1]>1e-5)printf("%lf\n", summess[ths][i][j][d1]);
					be[i][j][d1] = L_init[i][j][d1] + summess[ths][i][j][d1];
					
				}
			}
		}
	}
}

void init(Frame*frames)
{
	int height =frames[0].height , width =frames[0].width;
	cout << width << " " << height << endl;
	img_width = width;
	img_height = height;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			ux[i][j] = 0;
			for (int d = 0; d < k; d++) {
				L_init[i][j][d] = 0;
			}
		}
	}
	cv::Mat trans=cv::Mat(height, width, CV_8UC1);
	for (int i = 0; i <= 0; i++) {
		//cout << i << endl;
		Frame *frame = &frames[i];
		//cout << frame->cam->R<<"???\n";
		for (int d = 1; d <= k; d++) {
			for (int j = 1; j < 2; j++) {
				Frame *framet = &frames[j];
				//cout << framet->cam->R << "???\n";
				cv::Mat w_c0, w_c1, w_c2, w_c3;
				double x0, y0, x1, y1, x2, y2, x3, y3,d0, d1, d2, d3;
				frame->GetWorldCoordFrmImgCoord(0, 0, 1 / (d*delta), w_c0);
				frame->GetWorldCoordFrmImgCoord(height, 0, 1 / (d*delta), w_c1);
				frame->GetWorldCoordFrmImgCoord(0, width, 1 / (d*delta), w_c2);
				frame->GetWorldCoordFrmImgCoord(height, width, 1 / (d*delta), w_c3);
				framet->GetImgCoordFrmWorldCoord(x0, y0, d0, w_c0);
				framet->GetImgCoordFrmWorldCoord(x1, y1, d1, w_c1);
				framet->GetImgCoordFrmWorldCoord(x2, y2, d2, w_c2);
				framet->GetImgCoordFrmWorldCoord(x3, y3, d3, w_c3);
				cout << x0 << " " << y0 << " \n" << x1 << " " << y1 << "\n" << x2 << " " << y2 <<"\n"<< x3 << " " << y3 << endl;
				//cout << d0 << " " << d1 << " " << d2 << " " << d3 << endl;
				double ax = (x1*d1 - x0 * d0) / height, bx = (x2*d2 - x0 * d0) / width;
				double ay = (y1*d1 - y0 * d0) / height, by = (y2*d2 - y0 * d0) / width;
				double ad = (d1 - d0) / height, bd = (d2 - d0) / width;
				//cout << tmpd << endl;1
				//printf("----%lf %lf %lf %lf\n", ax, bx,ay,by);
				for (int x = 0; x < height; x++) {
					for (int y = 0; y < width; y++) {
						trans.at<uchar>(x, y) = 0;
					}
				}
				for (int x = 0; x < height; x++) {
					for (int y = 0; y <width; y++) {
						
						double tmp = 0;
						double td = (ad*x + bd * y + d0);
						//printf("%d %d %lf\n", x, y, td);
						int tx = (int)((ax*x + bx * y + x0*d0)/td+0.5);
						int ty = (int)((ay*x + by * y + y0*d0)/td+0.5);
						//if (x == 20 && y == 24) {
						//	printf("~~~%d %d\n", tx, ty);
						//}
						//printf("%d %d %d %d\n", x, y, tx, ty);
						if (tx >= 0 && tx < height&&ty >= 0 && ty < width) {
							trans.at<uchar>(tx, ty) = frame->gray_img.at<uchar>(x, y);
							tmp = abs(framet->gray_img.at<uchar>(tx, ty) - frame->gray_img.at<uchar>(x, y));
							tmp = 1 / (1 + tmp * tmp);
						}
						L_init[x][y][d - 1] += tmp;
						ux[x][y] = max(ux[x][y], L_init[x][y][d - 1]);
					}
				}
			}
			char file_name[20];
			sprintf(file_name, "%d.bmp", d);
			imwrite(file_name, trans);
		}
		//printf("%lf\n", ux);
		int *data = new int[width*height*k],*res= new int[width*height];
		int tot = 0,sum=0;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				for (int d = 0; d < k; d++) {
					if (ux[i][j] < 1e-8)
						L_init[i][j][d] = 0;
					else
						L_init[i][j][d] = 1 - 1.0 / ux[i][j] * L_init[i][j][d];
					data[tot++] = (int)(L_init[i][j][d] * 100);
					sum += data[tot - 1];
					//printf("!%d\n", data[tot - 1]);
				}
				double sumg = 0, sum = 0;
				if (i > 0) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i - 1, j));
				if (j > 0) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i, j-1));
				if (i < height-1) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i + 1, j));
				if (j < width-1) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i, j+1));
				u2[i][j] = sumg / sum;
			}
		}
		printf("%d %d %d\n",tot,sum,data[5174400]);
		gc(data,res);
		//bp(frame->gray_img);
		//cout << "bp done\n" << endl;
		tot = 0;
		for (int x = 0; x < height; x++) {
			for (int y = 0; y < width; y++) {
				//int ans = 0;
				//for (int d = 1; d < k; d++) {
				//	if (be[x][y][d] < be[x][y][ans]) {
				//		ans = d;
				//	}
				//}
				//if (ans != 0) printf("%d %d\n", x, y);
				frame->dep_img.at<uchar>(x, y)= res[tot++]*255/k;
			}
		}
		cv::imwrite("dep.bmp",frame->dep_img);
	}
}

void makepic()
{
	Camera * cam1 = new Camera(1, 1, 0, 0);
	Camera * cam2 = new Camera(1, 1, 0, 0);
	Frame *fm1 = new Frame();
	Frame *fm2 = new Frame();
	fm1->cam = cam1;
	fm2->cam = cam2;
	
	FILE* fp = fopen("test.txt", "r");

	cam1->Read(fp);
	cam2->Read(fp);

	cv::Mat pic(100, 100, CV_8UC1);
	cv::Mat rpic(100, 100, CV_8UC1);
	int n = 100, m = 100;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			pic.at<uchar>(i, j) = 240;
			rpic.at<uchar>(i, j) = 240;
		}
	}
	for (int i = 20; i < 40; i++) {
		for (int j = 20; j < 40; j++) {
			pic.at<uchar>(i, j) = 80;
		}
	}
	for (int i = 70; i < 90; i++) {
		for (int j = 70; j < 90; j++) {
			pic.at<uchar>(i, j) = 160;
		}
	}
	cv::imwrite("unrot.bmp", pic);
	fm1->gray_img = pic;
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			cv::Mat w_c;
			fm2->out = 0;
			fm1->GetWorldCoordFrmImgCoord(i, j, pic.at<uchar>(i, j),w_c);
			double a, b, c;
			fm2->GetImgCoordFrmWorldCoord(a, b, c, w_c);
			int x = a + (a>0?0.5:-0.5), y = b + (b > 0 ? 0.5 : -0.5), z = c + (c > 0 ? 0.5 : -0.5);
			if (i == 20 && j == 24) printf("???%d %d\n", x, y);
			if (x >= 0 && x < 100 && y>=0 && y < 100) {
				rpic.at<uchar>(x, y) = pic.at<uchar>(i, j);
			}
		}
	}
	cv::imwrite("rot.bmp", rpic);
}

int main()
{
	//makepic();
	//return 0;
	FILE* actfile = fopen("DATA/myact.txt", "r");
	int picnums;
	fscanf(actfile, "%d", &picnums);
	Frame* frames = new Frame[picnums];
	double fx,fy,cx,cy,s;
	fscanf(actfile,"%lf%lf%lf%lf%lf", &fx, &fy, &cx, &cy, &s);
	for (int i = 0; i < picnums; i++) {
		char filename[100];
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
	system("pause");
}