#include<iostream>
#include<core/core.hpp>
#include<highgui/highgui.hpp>  
#include"camera.h"
#include"frame.h"
#include "gcoptimization.h"

using namespace std;


int img_width, img_height;
#define maxd 30
double delta = 0.015/maxd;
int k = maxd;
double l_init[1000][1000][maxd];
double ux[1000][1000], u2[1000 * 1000];

double lum[1000 * 1000];

int src(int x, int y, int label)
{
	return (x * img_width + y)*k + label;
}

struct fordatafn {
	int numlab;
	int *data;
};


double sqr(double x)
{
	return x * x;
}

int datafn(int p, int l, void *data)
{
	
	fordatafn *mydata = (fordatafn *)data;
	int numlab = mydata->numlab;
	return(mydata->data[p*numlab + l]);
}

int smoothfn(int p1, int p2, int l1, int l2)
{
	double t = (u2[p1]+u2[p2])/ (sqr(lum[p1] - lum[p2]) + 1)*std::abs(l1 - l2);
	t *= 150;
	//if (t > 10) printf("%lf\n", t);
	return t;
}

void gc(int * data,int* &result)
{
	try {
		gcoptimizationgridgraph *gc = new gcoptimizationgridgraph(img_width, img_height, k);
		fordatafn tofn;
		printf("%d\n", data[5174400]);
		tofn.data = data;
		tofn.numlab = k;
		gc->setdatacost(&datafn,&tofn);
		gc->setsmoothcost(&smoothfn);
		printf("\nbefore optimization energy is %lld", gc->compute_energy());
		gc->swap(2);// run expansion for 2 iterations. for swap use gc->swap(num_iterations);
		printf("\nafter optimization energy is %lld", gc->compute_energy());

		for (int i = 0; i < img_width*img_height; i++)
			result[i] = gc->whatlabel(i);
		delete gc;
	}
	catch (gcexception e) {
		e.report();
		while (1);
	}
	/*double *nl = new double[img_width*img_height*k];
	double *x = l_init, *y = nl;
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
	l_init = x;
	delete y;
	int num_pixels = img_height * img_width;
	result = new int[num_pixels];   // stores result of optimization
	gcoptimizationgridgraph *gc = new gcoptimizationgridgraph(img_width, img_height, k);
	for (int i = 0; i < img_height*img_width; i++) {
		for (int j = 0; j < k; j++) {
			gc->setdatacost(i, j, l_init[i*k+j]);
		}
	}
	for (int l1 = 0; l1 < k; l1++)
		for (int l2 = 0; l2 < k; l2++) {
			int cost = (l1 - l2)*(l1 - l2);
			gc->setsmoothcost(l1, l2, 0);
		}
	cout << gc->compute_energy();
	gc->expansion(5);
	cout << gc->compute_energy();
	for (int i = 0; i < num_pixels; i++)
		result[i] = gc->whatlabel(i);*/
}


/*
double calsum(int x, int y, int x2, int y2, int iter, int d1, int d2,cv::mat &frame,int dir)
{
	int it = frame.at<uchar>(x, y) - frame.at<uchar>(x2, y2);
	double sum = l_init[x2][y2][d2] + abs(d1 - d2)*u2[x*img_width+y] / abs(it+1e-7);
	sum += summess[iter ^ 1][x2][y2][d2] - mess[iter ^ 1][x2][y2][d2][dir];
	return sum;
}

void bp(cv::mat &frame)
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
					be[i][j][d1] = l_init[i][j][d1] + summess[ths][i][j][d1];
					
				}
			}
		}
	}
}*/

void init(frame*frames)
{
	int height =frames[0].height , width =frames[0].width;
	cout << width << " " << height << endl;
	img_width = width;
	img_height = height;
	
	cv::mat trans=cv::mat(height, width, cv_8uc1);
	for (int i = 0; i < 140; i+=10) {
		for (int x = 0; x < height; x++) {
			for (int y = 0; y < width; y++) {
				ux[x][y] = 0;
				lum[x*width + y] = frames[i].gray_img.at<uchar>(x, y);
				for (int d = 0; d < k; d++) {
					l_init[x][y][d] = 0;
				}
			}
		}
		//cout << i << endl;
		frame *frame = &frames[i];
		//cout << frame->cam->r<<"???\n";
		for (int d = 1; d <= k; d++) {
			for (int j = max(i-30,0); j < min(i+30,140); j+=3) {
				frame *framet = &frames[j];
				uchar*pic1 = frame->gray_img.data;
				uchar*pic2 = framet->gray_img.data;
				//cout << framet->cam->r << "???\n";
				cv::mat w_c0, w_c1, w_c2, w_c3;
				double x0, y0, x1, y1, x2, y2, x3, y3,d0, d1, d2, d3;
				frame->getworldcoordfrmimgcoord(0, 0, 1 / (d*delta), w_c0);
				frame->getworldcoordfrmimgcoord(width, 0, 1 / (d*delta), w_c1);
				frame->getworldcoordfrmimgcoord(0, height, 1 / (d*delta), w_c2);
				frame->getworldcoordfrmimgcoord(width, height, 1 / (d*delta), w_c3);
				framet->getimgcoordfrmworldcoord(x0, y0, d0, w_c0);
				framet->getimgcoordfrmworldcoord(x1, y1, d1, w_c1);
				framet->getimgcoordfrmworldcoord(x2, y2, d2, w_c2);
				framet->getimgcoordfrmworldcoord(x3, y3, d3, w_c3);
				cout << x0 << " " << y0 << " \n" << x1 << " " << y1 << "\n" << x2 << " " << y2 <<"\n"<< x3 << " " << y3 << endl;
				//cout << d0 << " " << d1 << " " << d2 << " " << d3 << endl;
				double ax = (x1*d1 - x0 * d0) / width, bx = (x2*d2 - x0 * d0) / height;
				double ay = (y1*d1 - y0 * d0) / width, by = (y2*d2 - y0 * d0) / height;
				double ad = (d1 - d0) / width, bd = (d2 - d0) / height;
				//cout << tmpd << endl;1
				//printf("----%lf %lf %lf %lf\n", ax, bx,ay,by);
				/*for (int x = 0; x < height; x++) {
					for (int y = 0; y < width; y++) {
						trans.at<uchar>(x, y) = 0;
					}
				}*/
				for (int x = 0; x < height; x++) {
					for (int y = 0; y <width; y++) {
						
						double tmp = 0;
						double td = (ad*y + bd * x + d0);
						//printf("%d %d %lf\n", x, y, td);
						int ty = (int)((ax*y + bx * x + x0*d0)/td+0.5);
						int tx = (int)((ay*y + by * x + y0*d0)/td+0.5);
						//if (x == 20 && y == 24) {
						//	printf("~~~%d %d\n", tx, ty);
						//}
						//printf("%d %d %d %d\n", x, y, tx, ty);
						if (tx >= 0 && tx < height&&ty >= 0 && ty < width) {
							//trans.at<uchar>(tx, ty) = frame->gray_img.at<uchar>(x, y);
							//tmp = abs(framet->gray_img.at<uchar>(tx, ty) - frame->gray_img.at<uchar>(x, y));
							tmp = abs(pic1[x*width + y] - pic2[tx*width + ty]);
							tmp = 1 / (1 + tmp * tmp);
						}
						l_init[x][y][d - 1] += tmp;
						ux[x][y] = max(ux[x][y], l_init[x][y][d - 1]);
					}
				}
				
			}
		}
		//printf("%lf\n", ux);
		int *data = new int[width*height*k],*res= new int[width*height];
		int tot = 0,sum=0;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				for (int d = 0; d < k; d++) {
					if (ux[i][j] < 1e-8)
						l_init[i][j][d] = 0;
					else
						l_init[i][j][d] = 1 - 1.0 / ux[i][j] * l_init[i][j][d];
					data[tot++] = (int)(l_init[i][j][d] * 100);
					sum += data[tot - 1];
					//printf("!%d\n", data[tot - 1]);
				}
				double sumg = 0, sum = 0;
				if (i > 0) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i - 1, j)+1);
				if (j > 0) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i, j-1)+1);
				if (i < height-1) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i + 1, j)+1);
				if (j < width-1) sumg += 1, sum += 1/sqr(frame->gray_img.at<uchar>(i, j) - frame->gray_img.at<uchar>(i, j+1)+1);
				u2[i*width+j] = sumg / sum;
			}
		}
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
				frame->dep_img.at<uchar>(x, y)= res[tot++]*220/k+35;
			}
		}
		char file_name[20];
		sprintf(file_name, "dep%d.bmp",i);
		cv::imwrite(file_name,frame->dep_img);
	}
}

void makepic()
{
	camera * cam1 = new camera(1, 1, 0, 0);
	camera * cam2 = new camera(1, 1, 0, 0);
	frame *fm1 = new frame();
	frame *fm2 = new frame();	
	fm1->cam = cam1;
	fm2->cam = cam2;
	
	file* fp = fopen("test.txt", "r");

	cam1->read(fp);
	cam2->read(fp);

	cv::mat pic(100, 100, cv_8uc1);
	cv::mat rpic(100, 100, cv_8uc1);
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
			cv::mat w_c;
			fm2->out = 0;
			fm1->getworldcoordfrmimgcoord(i, j, pic.at<uchar>(i, j),w_c);
			double a, b, c;
			fm2->getimgcoordfrmworldcoord(a, b, c, w_c);
			int x = a + (a>0?0.5:-0.5), y = b + (b > 0 ? 0.5 : -0.5), z = c + (c > 0 ? 0.5 : -0.5);
			if (i == 20 && j == 24) printf("???%d %d\n", x, y);
			if (x >= 0 && x < 100 && y>=0 && y < 100) {
				rpic.at<uchar>(x, y) = pic.at<uchar>(i, j);
			}
		}
	}
	cv::imwrite("rot.bmp", rpic);
}

void depthrecover(frame*frames)
{
	init(frames);
}

void calcpoint(frame*frames)
{
	frame*frame = &frames[70];
	frame->dep_img = cv::imread("dep70.bmp", cv_8uc1);
	int height = frame->height, width = frame->width;
	file*out=fopen("points.txt", "w");
	for (int i = 0; i < height; i+=2) {
		for (int j = 0; j < width; j+=2) {
			if (frame->dep_img.data[i*width + j] < 50) continue;
			cv::mat w_c;
			frame->getworldcoordfrmimgcoord(j, i, 1/(k*((frame->dep_img.data[i*width + j]-35)/220.0+1)*delta), w_c);
			//fprintf(out,"%d %d %d\n", i, j, frame->dep_img.data[i*width + j]);
			fprintf(out, "%d %d %d\n", (int)(w_c.at<float>(0,0)+0.5), (int)(w_c.at<float>(1, 0) + 0.5), (int)(w_c.at<float>(2, 0) + 0.5));
		}
	}
	fclose(out);
}

int main()
{
	file* actfile = fopen("data/myact.txt", "r");
	int picnums;
	fscanf(actfile, "%d", &picnums);
	frame* frames = new frame[picnums];
	double fx, fy, cx, cy, s;
	fscanf(actfile, "%lf%lf%lf%lf%lf", &fx, &fy, &cx, &cy, &s);
	for (int i = 0; i < picnums; i++) {
		char filename[100];
		fscanf(actfile, "%s", filename);
		camera*cam = new camera(fx, fy, cx, cy);
		fscanf(actfile, "%lf", &s);
		cam->read(actfile);
		frames[i].cam = cam;
		sprintf(filename, "data/test%04d.jpg", i);
		frames[i].readimg(filename);
		fscanf(actfile, "%lf,%lf,%lf,%lf", &s, &s, &s, &s);
		fscanf(actfile, "%s", filename);
	}
	//depthrecover(frames);
	calcpoint(frames);
	system("pause");
}
