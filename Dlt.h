#pragma once
#ifndef DLT
#define DLT
#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

class Dlt
{
public:
	Dlt();
	Dlt(double * in_x, double * in_y, double * in_xt, double * in_yt);
	~Dlt();

	Mat getMatH();
	void normalize(Mat& o_T, Mat& t_T, Dlt& normal);
	double getValue(char c, bool trans, int index);
	
private:
	double * x;
	double * y;
	double * xt;
	double * yt;

	Mat mat_Init();
};


#endif // !DLT
