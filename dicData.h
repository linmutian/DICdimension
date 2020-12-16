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

enum Progression { linear, powerLaw }; // ����/����
const double pi = 3.14159265358979323846;
// ���(x-A)^2+(y-B)^2=R^2
bool circFit(vector<double>x, vector<double>y, double &A, double &B, double &R);
// ��� y=kx+b
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
	/// ���ú���
	void cartesian2polar(); // x,y,z->xs,y
	void calcRange(); //
	void filter();
	void calcDf();
	bool readConfig(string fileName);
	bool writeData(string fileName);

	/// ���ò��� 
	double rangeFac = 0.02;  // ��Χϵ��
	double sigThr = 1.00;	// Ӧ����ֵ
	int nNum = 30; // ������ϵ����ݶ���
	int nBegin = 2; // ��������������ÿ��ά��n�ȷ֣�n��2��ʼ
	Progression pro = linear;

	/// ������
	vector<double> xVec;
	vector<double> yVec;
	vector<double> zVec; 
	vector<double> sigVec; // ��Ӧ��
	vector<double> xsVec;
	// ɸѡ��
	vector<double> xsFilVec;
	vector<double> yFilVec;
	// �������
	vector<double> epsVec; // ��
	vector<double> NVec; // N(��)
	vector<double> lnEpsVec;
	vector<double> lnNVec;
	// ���㷶Χ
	double xsRangeL;
	double xsRangeH;
	double yRangeL;
	double yRangeH;
};