#include "dicData.h"

bool circFit(vector<double> x, vector<double> y, double & cenX, double & cenY, double & R)
{
	if (!x.size() || x.size() != y.size())
		return false;
	double x1 = 0.0;
	double x2 = 0.0;
	double x3 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;
	double y3 = 0.0;
	double x1y1 = 0.0;
	double x1y2 = 0.0;
	double x2y1 = 0.0;
	size_t N = x.size();
	for (size_t i = 0; i < N; i++)
	{
		x1 += x[i];
		x2 += x[i] * x[i];
		x3 += x[i] * x[i] * x[i];
		y1 += y[i];
		y2 += y[i] * y[i];
		y3 += y[i] * y[i] * y[i];
		x1y1 += x[i] * y[i];
		x1y2 += x[i] * y[i] * y[i];
		x2y1 += x[i] * x[i] * y[i];
	}
	double C = N * x2 - x1 * x1;
	double D = N * x1y1 - x1 * y1;
	double E = N * x3 + N * x1y2 - (x2 + y2) * x1;
	double G = N * y2 - y1 * y1;
	double H = N * x2y1 + N * y3 - (x2 + y2) * y1;
	double a = (H * D - E * G) / (C * G - D * D);
	double b = (H * C - E * D) / (D * D - G * C);
	double c = -(a * x1 + b * y1 + x2 + y2) / N;
	cenX = a / (-2);
	cenY = b / (-2);
	R = sqrt(a * a + b * b - 4 * c) / 2;
	return true;
}

double lineFit(vector<double> x, vector<double> y, double & k, double & b)
{
	double av_x, av_y; 
	double L_xx, L_yy, L_xy;
	int num = x.size();
	av_x = 0; //X的平均值
	av_y = 0; //Y的平均值
	L_xx = 0; //Lxx
	L_yy = 0; //Lyy
	L_xy = 0; //Lxy

	for (int i = 0; i < num; i++) 
	{
		av_x += x[i];
		av_y += y[i];
	}
	av_x /= num;
	av_y /= num;

	for (int i = 0; i < num; i++) 
	{
		L_xx += (x[i] - av_x)*(x[i] - av_x);
		L_yy += (y[i] - av_y)*(y[i] - av_y);
		L_xy += (x[i] - av_x)*(y[i] - av_y);
	}
	k = L_xy / L_xx;
	b = av_y - L_xy * av_x / L_xx;
	double r = L_xy / sqrt(L_xx*L_yy);
	double r2 = pow(r, 2);
	return r2;
}

void init(vector<vector<bool>>& boxes, int n)
{
	vector<bool> row(n, false);
	boxes.clear();
	for (size_t i = 0; i < n; i++)
	{
		boxes.push_back(row);
	}
}

dicData::dicData(string path)
{
	ifstream fread(path);
	if (fread.fail())
	{
		cout << "[error-0] 数据文件打开失败！ " << endl;
		return;
	}
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double sig = 0.0;
	while (!fread.eof())
	{
		fread >> x >> y >> z >> sig;
		xVec.push_back(x);
		yVec.push_back(y);
		zVec.push_back(z);
		sigVec.push_back(sig);
	}
}

void dicData::set_nBegin(int nb) { nBegin = nb; }

void dicData::set_nNum(int nn) { nNum = nn; }

void dicData::setMode2Line() { pro = linear; }

void dicData::setMode2Power() { pro = powerLaw; }

void dicData::run()
{
	cartesian2polar();
	filter();
	calcRange();
	calcDf();
}

void dicData::runWithConfig()
{
	if (readConfig("config.txt"))
	{
		cartesian2polar();
		filter();
		calcRange();
		calcDf();
		if(writeData("dimensionData.csv"))
			cout << "[info] 拟合数据已写入dimensionData.csv中" << endl;
	}
	return;
}

void dicData::cartesian2polar()
{
	// 拟合圆
	double cenX = 0.0;
	double cenY = 0.0;
	double R = 0.0;
	circFit(xVec, zVec, cenX, cenY, R);
#ifdef DEBUG
	cout << "[debug] 拟合结果：(x-" << cenX << ")^2+(y-" << cenY << ")^2 = " << R << "^2" << endl;
#endif // DEBUG

	//转换为极坐标
	vector<double> rVec;
	vector<double> phiVec; // φ
	for (size_t i = 0; i < xVec.size(); i++)
	{
		rVec.push_back(sqrt(pow(xVec[i] - cenX, 2) + pow(zVec[i] - cenY, 2)));
		phiVec.push_back(atan2(xVec[i] - cenX, zVec[i] - cenY));
	}

	// φ调整至连续 
	double minDistSum = DBL_MAX; // 最小距离和
	double phiMid = 0.0; // 估计的数据点中心
	for (double phiM = -pi; phiM < pi; phiM += 0.1*pi)
	{
		double distSum = 0.0;
		for (double &phi : phiVec)
		{
			if (phi > phiM + pi)
				distSum += abs(phi - pi - phiM);
			else if (phi < phiM - pi)
				distSum += abs(phi + pi - phiM);
			else
				distSum += abs(phi - phiM);
		}
		if (distSum < minDistSum)
		{
			minDistSum = distSum;
			phiMid = phiM;
		}
	}
	for (double &phi : phiVec)
	{
		phi -= phiMid;
		if (phi >  pi)
			phi -= pi;
		else if (phi < - pi)
			phi += pi;
		xsVec.push_back(phi*R);
	}
}

void dicData::calcRange()
{
	xsRangeL = *std::min_element(xsVec.begin(), xsVec.end());
	xsRangeH = *std::max_element(xsVec.begin(), xsVec.end());
	yRangeL = *std::min_element(yVec.begin(), yVec.end());
	yRangeH = *std::max_element(yVec.begin(), yVec.end());

	double xDiff = xsRangeH - xsRangeL;
	double yDiff = yRangeH - yRangeL;
	xsRangeL -= xDiff * rangeFac;
	xsRangeH += xDiff * rangeFac;
	yRangeL -= yDiff * rangeFac;
	yRangeH += yDiff * rangeFac;
}

void dicData::filter()
{
	xsFilVec.clear();
	yFilVec.clear();
	for (size_t i = 0; i < xVec.size(); i++)
	{
		if (sigVec[i] >= sigThr)
		{
			xsFilVec.push_back(xsVec[i]);
			yFilVec.push_back(yVec[i]);
		}
	}
#ifdef DEBUG
	cout <<"[debug] " <<xsVec.size() << " filtered to " << xsFilVec.size() << endl;
#endif // DEBUG
}

void dicData::calcDf()
{
	vector<vector<bool>> boxes; // true - 盒子中包含此形状
	int N_eps = 0;
	int n = nBegin;
	for (int ii=0; ii < nNum; ii++) // 每个维度n等分
	{
		N_eps = 0;
		init(boxes, n);
		for (size_t i = 0; i < xsFilVec.size(); i++)
		{
			
			int row = floor((xsFilVec[i] - xsRangeL) / (xsRangeH - xsRangeL)*n);
			int col = floor((yFilVec[i] - yRangeL) / (yRangeH - yRangeL)*n);
			// (xsFilVec[i] - xsRangeL) / (xsRangeH - xsRangeL)应该在0~1之间，不可能越界
			boxes[row][col] = true;
		}
		for (auto row : boxes)
		{
			for (bool ele : row)
			{
				if (ele)
					N_eps++;
			}
		}
		epsVec.push_back(1.0/n);
		NVec.push_back(N_eps);
		lnEpsVec.push_back(-log(n));
		lnNVec.push_back(log(N_eps));

		if (pro == linear)
			n++;
		else
			n *= 2;
	}
	double k, b;
	double r2 = lineFit(lnEpsVec, lnNVec, k, b);
	cout << "[info] 拟合方程为：ln(N)= " << k << " × ln(ε) + " << b << "   , r^2 = " << r2 << endl;
	cout << "       即，分形维数 Df = " << -k << endl;
}

bool dicData::readConfig(string fileName)
{
	ifstream configIn(fileName);
	if (configIn)
	{
		int mode;
		configIn >> sigThr>> rangeFac>> nBegin >> nNum >> mode;
		if (mode == 2)
			setMode2Power();
		else
			setMode2Line();
		// 输出设置
		cout << "[info] 计算设置：σ_thr = " << sigThr << " , fr = " << rangeFac 
			<< " , n0 = " << nBegin << ", num_n = " << nNum << endl;
		cout << "[info] 序列模式：" << (mode == 2 ? "幂律" : "线性") << endl;

		configIn.close();
		return true;
	}
	else
	{
		cout << "[error-1] 配置文件无法读取！本exe同路径下应存在config.txt ！" << endl;
		return false;
	}
}

bool dicData::writeData(string fileName)
{
	ofstream writeDataStream(fileName);
	if (writeDataStream.fail())
	{
		cout << "[error-2] 输出文件无法打开!请检查其是否被占用！" << endl;
		return false;
	}
	writeDataStream << "ε,N,ln(ε),ln(N)" << endl;
	for (int i = 0; i < lnNVec.size(); i++)
	{
		writeDataStream << epsVec[i] << "," << NVec[i] << "," << lnEpsVec[i] << "," << lnNVec[i] << endl;
	}
	writeDataStream.close();
	return true;
}
