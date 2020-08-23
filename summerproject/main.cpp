#include<iostream>
#include<core/core.hpp>
#include<highgui/highgui.hpp>  
#include"camera.h"
#include"frame.h"
#include "GcOptimization.h"

using namespace std;


int img_width, img_height;
#define maxd 50
#define disp 0.015
double delta = disp / maxd;
int k = maxd;
double l_init[1000][1000][maxd];
double ux[1000][1000], u2[1000 * 1000];

double lum[1000 * 1000];
double eps = 50,Ws=10/disp;
Frame*nowFrame;

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
	double t = (u2[p1] + u2[p2]) / (sqr(lum[p1] - lum[p2]) + 1)*std::abs(l1 - l2);
	t *= 150;
	//if (t > 10) printf("%lf\n", t);
	return t;
}

void gc(int * data, int* &result)
{
	try {
		GCoptimizationGridGraph *gc = new GCoptimizationGridGraph(img_width, img_height, k);
		fordatafn tofn;
		printf("%d\n", data[5174400]);
		tofn.data = data;
		tofn.numlab = k;
		gc->setDataCost(&datafn, &tofn);
		gc->setSmoothCost(&smoothfn);
		printf("\nbefore optimization energy is %lld", gc->compute_energy());
		gc->swap(2);// run expansion for 2 iterations. for swap use gc->swap(num_iterations);
		printf("\nafter optimization energy is %lld", gc->compute_energy());

		for (int i = 0; i < img_width*img_height; i++)
			result[i] = gc->whatLabel(i);
		delete gc;
	}
	catch (GCException e) {
		e.Report();
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

float calcdis(int x1, int y1, int x2, int y2)
{
	int p1 = x1 * nowFrame->width + y1;
	int p2 = x2 * nowFrame->width + y2;
	int flag = 0;
	if (flag) {
		uchar*data = nowFrame->gray_img.data;
		return sqr(data[p1] - data[p2]);
	}
	else {
		uchar*data = nowFrame->color_img.data;
		return sqr(data[p1 * 3] - data[p2 * 3]) + sqr(data[p1 * 3 + 1] - data[p2 * 3 + 1]) + sqr(data[p1 * 3 + 2] - data[p2 * 3 + 2]);
	}
}

namespace bp {
	const int LEVELS = 5, ITER = 20;
	const float disthreshold = 0.05*disp;
	struct DATA {
		int width, height;
		float*dat;
		void set(int x, int y, int d,float v) {
			dat[(x*height + y)*maxd + d]=v;
		}
		void add(int x, int y, int d, float v) {
			dat[(x*height + y)*maxd + d] += v;
		}
		float *get(int x, int y, int d) {
			return &dat[(x*height + y)*maxd + d];
		}
		float *get(int x, int y) {
			return &dat[(x*height + y)*maxd];
		}
		void init(int _width,int _height) {
			width = _width;
			height = _height;
			dat = new float[width*height*maxd]();
		}
		void del() {
			delete[]dat;
		}
	}data[LEVELS],u[LEVELS],d[LEVELS],l[LEVELS],r[LEVELS];
	void dt(float *f,float lbd) {
		for (int q = 1; q < maxd; q++) {
			float pre = f[q - 1] + lbd;
			if (pre < f[q]) {
				f[q] = pre;
			}
		}
		for (int q = maxd - 2; q >= 0; q--) {
			float pre = f[q + 1] + lbd;
			if (pre < f[q]) {
				f[q] = pre;
			}
		}
	}
	void msg(float*s1, float*s2, float*s3, float*s4, float*dst,float lbd) {
		float mn = 1e9;
		for (int d = 0; d < maxd; d++) {
			dst[d] = s1[d] + s2[d] + s3[d] + s4[d];
			if (dst[d] < mn) {
				mn = dst[d];
			}
		}
		dt(dst, lbd);
		mn += disthreshold*lbd;
		for (int d = 0; d < maxd; d++) {
			if (mn < dst[d]) {
				dst[d] = mn;
			}
		}
		float val = 0;
		for (int d = 0; d < maxd; d++) {
			val += dst[d];
		}
		val /= maxd;
		for (int d = 0; d < maxd; d++) {
			dst[d] -= val;
		}
	}
	float calclambda(int x,int y,int x1,int y1) {
		int p1 = x * img_width + y, p2 = x1 * img_width + y1;
		return Ws*u2[p1]/(calcdis(x,y,x1,y1)+eps);
	}
	void BP(DATA&u, DATA&d, DATA&l, DATA&r, DATA&data,int ITER) {
		int width = data.width, height = data.height;
		for (int t = 0; t < ITER; t++) {
			for (int y = 1; y < height - 1; y++) {
				for (int x = ((y + t) & 1) + 1; x < width - 1; x += 2) {
					msg(u.get(x, y + 1), l.get(x + 1, y), r.get(x - 1, y), data.get(x, y), u.get(x, y), calclambda(x, y, x, y - 1));
					msg(d.get(x, y - 1), l.get(x + 1, y), r.get(x - 1, y), data.get(x, y), d.get(x, y), calclambda(x, y, x, y + 1));
					msg(u.get(x, y + 1), d.get(x, y - 1), r.get(x - 1, y), data.get(x, y), r.get(x, y), calclambda(x + 1, y, x, y));
					msg(u.get(x, y + 1), d.get(x, y - 1), l.get(x + 1, y), data.get(x, y), l.get(x, y), calclambda(x - 1, y, x, y));
				}
			}
		}

	}
	uchar* output(DATA&u, DATA&d, DATA&l, DATA&r, DATA&data){
		int width = data.width, height = data.height;
		uchar *dep = new uchar[width*height]();
		for (int y = 1; y < height - 1; y++) {
			for (int x = 1; x < width - 1; x++) {
				int best = 0;
				float best_val = 1e9;
				for (int de = 0; de < maxd; de++) {
					float val = *u.get(x, y + 1, de) +
						*d.get(x, y - 1, de) +
						*l.get(x + 1, y, de) +
						*r.get(x - 1, y, de) + 
						*data.get(x, y, de);
					if (val < best_val) {
						best_val = val;
						best = de;
					}
				}
				dep[x*height + y] = 1.0*best / maxd * 255;
			}
		}
		return dep;
	}
	void runbp(int width,int height,int id) {
		cout << "get into bp\n";
		data[0].init(width,height);
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int d = 0; d < maxd; d++) {
					data[0].set(i, j, d, l_init[i][j][d]);
				}
			}
		}
		for (int i = 1; i < LEVELS; i++) {
			data[i].init(data[i-1].width/2.0+0.6,data[i-1].height/2.0+0.6);
			cout << "~" << data[i - 1].width / 2 << " " << data[i - 1].height / 2 << endl;
			for (int x = 0; x < data[i-1].width; x++) {
				for (int y = 0; y < data[i-1].height; y++) {
					for (int d = 0; d < maxd; d++) {
						//cout << "!" << x << " " << y << endl;
						data[i].add(x/2, y/2, d, *data[i - 1].get(x, y, d));
					}
				}
			}
		}
		for (int i = LEVELS - 1; i >= 0; i--) {
			u[i].init(data[i].width, data[i].height);
			d[i].init(data[i].width, data[i].height);
			l[i].init(data[i].width, data[i].height);
			r[i].init(data[i].width, data[i].height);
			int width = data[i].width, height = data[i].height;
			if (i != LEVELS - 1) {
				for (int x = 0; x < width; x++) {
					for (int y = 0; y < height; y++) {
						for (int de = 0; de < maxd; de++) {
							u[i].set(x, y, de, *u[i + 1].get(x / 2, y / 2, de));
							d[i].set(x, y, de, *d[i + 1].get(x / 2, y / 2, de));
							l[i].set(x, y, de, *l[i + 1].get(x / 2, y / 2, de));
							r[i].set(x, y, de, *r[i + 1].get(x / 2, y / 2, de));
						}
					}
				}
			}
			u[i + 1].del();
			d[i + 1].del();
			l[i + 1].del();
			r[i + 1].del();
			BP(u[i], d[i], l[i], r[i], data[i], ITER);
		}
		cv::Mat depth=cv::Mat(width,height,CV_8UC1,output(u[0], d[0], l[0], r[0], data[0]));
		char file_name[20];
		sprintf(file_name, "depth%d.bmp", id);
		cv::imwrite(file_name, depth);
		cout << "bp end\n";
	}
}

void init(Frame*frames)
{
	int height = frames[0].height, width = frames[0].width;
	cout << width << " " << height << endl;
	img_width = width;
	img_height = height;

	cv::Mat trans = cv::Mat(height, width, CV_8UC1);
	for (int i = 20; i < 22; i ++) {
		printf("%d\n", i);
		for (int x = 0; x < height; x++) {
			for (int y = 0; y < width; y++) {
				ux[x][y] = 0;
				lum[x*width + y] = frames[i].gray_img.at<uchar>(x, y);
				for (int d = 0; d < k; d++) {
					l_init[x][y][d] = 0;
				}
			}
		}

		Frame *frame = &frames[i];
		nowFrame = &frames[i];
		for (int d = 1; d <= k; d++) {
			for (int j = i-15; j < i+15; j++) {
				Frame *framet = &frames[j];
				uchar*pic1 = frame->gray_img.data;
				uchar*pic2 = framet->gray_img.data;

				cv::Mat w_c0, w_c1, w_c2, w_c3;
				double x0, y0, x1, y1, x2, y2, x3, y3, d0, d1, d2, d3;
				frame->GetWorldCoordFrmImgCoord(0, 0, 1 / (d*delta), w_c0);
				frame->GetWorldCoordFrmImgCoord(width, 0, 1 / (d*delta), w_c1);
				frame->GetWorldCoordFrmImgCoord(0, height, 1 / (d*delta), w_c2);
				frame->GetWorldCoordFrmImgCoord(width, height, 1 / (d*delta), w_c3);
				framet->GetImgCoordFrmWorldCoord(x0, y0, d0, w_c0);
				framet->GetImgCoordFrmWorldCoord(x1, y1, d1, w_c1);
				framet->GetImgCoordFrmWorldCoord(x2, y2, d2, w_c2);
				framet->GetImgCoordFrmWorldCoord(x3, y3, d3, w_c3);
				double ax = (x1*d1 - x0 * d0) / width, bx = (x2*d2 - x0 * d0) / height;
				double ay = (y1*d1 - y0 * d0) / width, by = (y2*d2 - y0 * d0) / height;
				double ad = (d1 - d0) / width, bd = (d2 - d0) / height;

				for (int x = 0; x < height; x++) {
					for (int y = 0; y < width; y++) {

						double tmp = 0;
						double td = (ad*y + bd * x + d0);
						int ty = (int)((ax*y + bx * x + x0 * d0) / td + 0.5);
						int tx = (int)((ay*y + by * x + y0 * d0) / td + 0.5);

						if (tx >= 0 && tx < height&&ty >= 0 && ty < width) {
							tmp = abs(pic1[x*width + y] - pic2[tx*width + ty]);
							tmp = 10 / (10 + tmp * tmp);
						}
						l_init[x][y][d - 1] += tmp;
						ux[x][y] = max(ux[x][y], l_init[x][y][d - 1]);
					}
				}

			}
		}

		int *data = new int[width*height*k], *res = new int[width*height];
		int tot = 0, sum = 0;
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				for (int d = 0; d < k; d++) {
					if (ux[i][j] < 1e-8)
						l_init[i][j][d] = 1;
					else
						l_init[i][j][d] = 1 - 1.0 / ux[i][j] * l_init[i][j][d];
					data[tot++] = (int)(l_init[i][j][d] * 100);
					sum += data[tot - 1];
				}
				double sumg = 0, sum = 0;
				if (i > 0) sumg += 1, sum += 1 / (calcdis(i, j, i - 1, j) + eps);
				if (j > 0) sumg += 1, sum += 1 / (calcdis(i, j, i, j - 1) + eps);
				if (i < height - 1) sumg += 1, sum += 1 / (calcdis(i, j, i + 1, j) +eps);
				if (j < width - 1) sumg += 1, sum += 1 / (calcdis(i, j, i, j + 1) +eps);
				u2[i*width + j] = sumg / sum;
			}
		}
		/*gc(data, res);
		tot = 0;
		for (int x = 0; x < height; x++) {
			for (int y = 0; y < width; y++) {
				frame->dep_img.at<uchar>(x, y) = res[tot++] * 220 / k + 35;
			}
		}*/
		/*uchar*predepth = new uchar[height*width];
		for (int x = 0; x < height; x++) {
			for (int y = 0; y < width; y++) {
				int ansd = 0; float ans = 1e9;
				for (int d = 0; d < maxd; d++) {
					if (l_init[x][y][d] < ans) {
						ans = l_init[x][y][d];
						ansd = d;
					}
				}
				predepth[x*width + y] = 1.0*ansd/maxd*255;
			}
		}
		char file_name[20];
		sprintf(file_name, "predep%d.bmp", i);
		cv::Mat predep = cv::Mat(height, width, CV_8UC1, predepth);
		cv::imwrite(file_name, predep);*/
		bp::runbp(img_height,img_width,i);
		
	}
}

void depthrecover(Frame*frames)
{
	init(frames);
}

void calcpoint(Frame*frames)
{
	Frame*frame = &frames[70];
	frame->dep_img = cv::imread("dep70.bmp", CV_8UC1);
	int height = frame->height, width = frame->width;
	FILE*out = fopen("points.txt", "w");
	for (int i = 0; i < height; i += 2) {
		for (int j = 0; j < width; j += 2) {
			if (frame->dep_img.data[i*width + j] < 50) continue;
			cv::Mat w_c;
			frame->GetWorldCoordFrmImgCoord(j, i, 1 / (k*((frame->dep_img.data[i*width + j] - 35) / 220.0 + 1)*delta), w_c);
			//fprintf(out,"%d %d %d\n", i, j, frame->dep_img.data[i*width + j]);
			fprintf(out, "%d %d %d\n", (int)(w_c.at<float>(0, 0) + 0.5), (int)(w_c.at<float>(1, 0) + 0.5), (int)(w_c.at<float>(2, 0) + 0.5));
		}
	}
	fclose(out);
}

int main()
{
	FILE* actfile = fopen("data/myact.txt", "r");
	int picnums;
	fscanf(actfile, "%d", &picnums);
	Frame* frames = new Frame[picnums];
	double fx, fy, cx, cy, s;
	fscanf(actfile, "%lf%lf%lf%lf%lf", &fx, &fy, &cx, &cy, &s);
	for (int i = 0; i < picnums; i++) {
		char filename[100];
		fscanf(actfile, "%s", filename);
		Camera*cam = new Camera(fx, fy, cx, cy);
		fscanf(actfile, "%lf", &s);
		cam->Read(actfile);
		frames[i].cam = cam;
		sprintf(filename, "data/test%04d.jpg", i);
		frames[i].readimg(filename);
		fscanf(actfile, "%lf,%lf,%lf,%lf", &s, &s, &s, &s);
		fscanf(actfile, "%s", filename);
	}
	depthrecover(frames);
	//calcdepth(frame);
	//calcpoint(frames);
	system("pause");
}