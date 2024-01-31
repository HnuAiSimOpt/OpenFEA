#pragma once
#include<string>
#include<vector>
#include<sstream>
#include<math.h>

#include <cstring>
#include <cstdarg>
#include <map>
#include <set>
using namespace std;

namespace CAE
{
	class DataConvertor
	{ 
	private:
		const string* m_dataPtr = NULL;		/*指向待转换的数据行*/
		/*数据声明为常量指针，在转换器内不能对数据行进行修改*/

		int m_transPtr = 0;		/*当前转换位置下标*/
	public:
		DataConvertor();
		void	inputDataLine(string& dataline);
		void	ignore(int length);
		int		convertToInt(string& target, int length);
		int		convertToInt(int length);
		double  convertToDouble(int length);
		double  convertToDouble(string& target, int length);
		string  convertToString(int length);
		string	convertToString(int a, int divideSymbsASCII);

	private:
		void __check(int& length);
		int stringToInt(string& target);
		double stringToDouble(string& target);
	};
}