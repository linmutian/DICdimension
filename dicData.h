#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
//#define DEBUG

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

enum Progression { linear, powerLaw }; // 线性/幂律
const double pi = 3.14159265358979323846;
// 拟合(x-A)^2+(y-B)^2=R^2
bool circFit(vector<double>x, vector<double>y, double &A, double &B, double &R);
// 拟合 y=kx+b
double lineFit(vector<double>x, vector<double>y, double &k, double &b);
void init(vector<vector<bool>>& boxes, int n);

class dicData
{
public:
	dicData(string path);
	void set_nBegin(int nb);
	void set_nNum(int nn);
	void setMode2Line();
	void setMode2Power();
	void run();
	void runWithConfig();

private:
	/// 调用函数
	void cartesian2polar(); // x,y,z->xs,y
	void calcRange(); //
	void filter();
	void calcDf();
	bool readConfig(string fileName);
	bool writeData(string fileName);

	/// 设置参数 
	double rangeFac = 0.02;  // 范围系数
	double sigThr = 1.00;	// 应力阈值
	int nNum = 30; // 用于拟合的数据对数
	int nBegin = 2; // 将整个数据区域每个维度n等分，n从2开始
	Progression pro = linear;

	/// 数据区
	vector<double> xVec;
	vector<double> yVec;
	vector<double> zVec; 
	vector<double> sigVec; // σ应力
	vector<double> xsVec;
	// 筛选后
	vector<double> xsFilVec;
	vector<double> yFilVec;
	// 拟合数据
	vector<double> epsVec; // ε
	vector<double> NVec; // N(ε)
	vector<double> lnEpsVec;
	vector<double> lnNVec;
	// 计算范围
	double xsRangeL;
	double xsRangeH;
	double yRangeL;
	double yRangeH;
};