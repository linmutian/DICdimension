# DICdimension
用于计算特定的DIC数据的局部化带分形维数

一、输入与设置
1.1 数据要求
本程序用于计算DIC数据的分形维数。程序输入为一个txt文档，应有四列数据，分别为空间坐标x，y，z和一个应变值。其中应变将用于筛选局部化带内的点。
其中，点(x,y,z)是圆柱试件的侧面点。圆柱体中心轴应近似与y轴平行。
 
1.2 设置文件
设置文件为与程序处于同一路径下的config.txt文件。文件开头为5个以空格的数字，如1.0 0.02 4 30 1，后面的文字无任何影响。
(1) 第一个数字为$σ_thr$ ，即用于判断是否位于局部化带内的应力阈值；   
(2) 第二个数字为边框宽度系数$f_r$，见2.2；
(3) 第三到五个数字代表n的取值方法。第三个数字为n0，通常取2；第四个数字为$num_n$；第五个数字为2时按幂律取，否则按线性取，见2.3。

![](http://latex.codecogs.com/svg.latex?$σ_thr$)
