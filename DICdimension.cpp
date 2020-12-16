#include"dicData.h"
int main(int argc, char *argv[])
{
	if (argc > 1)
	{
		string dataFileName = argv[1];
		dicData d(dataFileName);
		d.runWithConfig();
	}
	else
	{
		cout << "[info] ����ģʽ" << endl;
		dicData d("D:\\work\\dic2.txt");
		d.run();
	}
	std::cin.get();
	return 0;
}