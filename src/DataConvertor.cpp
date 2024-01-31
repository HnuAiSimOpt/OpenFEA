#include "include/DataConvertor.h"
namespace CAE
{
	DataConvertor::DataConvertor()
    {}
    
	/* *************************Function*********************************
	函数名称：inputDataLine
	函数功能：接受需要转换的一行数据，并将当前转换位置指针移动到行首部。
	输入参数：待转换的数据行
	作者：GDZ
	日期：2021年7月14日 */
    void DataConvertor::inputDataLine(string &dataline)

	{
		m_dataPtr = &dataline;	//指向新的待转换数据
		m_transPtr = 0;	//重置转换位置
	}
	/* *************************Function*********************************
	函数名称：ignore
	函数功能：将从当前转换位置开始，长度为length的字符忽略
	输入参数：“擦除"字符的长度length
	作者：GDZ
	日期：2021年7月14日 */
	void DataConvertor::ignore(int length)
	{
		__check(length);	/* 确保长度不越界 */
		m_transPtr += length;
	}

	/* *************************Function*********************************
	函数名称：convertToInt
	函数功能：将从当前转换位置开始，且长度为length的字符子串转换为Int型数据，并将
			  转换位置指针向后移动length个位置，后续转换时将从新的转换位置开始。
	输入参数：转换字符串的长度length
	作者：GDZ
	日期：2021年7月14日 */
	int DataConvertor::convertToInt(int length)
	{
		__check(length);	/* 确保长度不越界 */
		/* 获取待转换子串 */
		string convertData = (*m_dataPtr).substr(m_transPtr, length);
		m_transPtr += length;			/* 移动转换位置，用于下次转换 */
		return stringToInt(convertData);	/* 转换子串为Int型并返回 */
	}

	int DataConvertor::convertToInt(string& target, int length)
	{
		length = length < target.size() ? length : target.size();	/* 确保长度不越界 */
		string convertData = target.substr(0, length);/* 获取待转换子串 */
		return stringToInt(convertData);	/* 转换子串为Int型并返回 */
	}
	/* *************************Function*******************************
	函数名称：convertToDouble
	函数功能：将从当前转换位置开始，且长度为length的字符子串转换为Double型数据，
			并将转换位置指针向后移动length个位置，后续转换时将从新的转换位置开始。
	输入参数：转换字符串的长度length
	作者：GDZ
	日期：2021年7月14日 */
	double DataConvertor::convertToDouble(int length)
	{
		__check(length);
		/* 获取待转换子串 */
		string convertData = (*m_dataPtr).substr(m_transPtr, length);
		/* 移动转换位置，用于下次转换 */
		m_transPtr += length;
		/* 转换子串为Double型并返回 */
		return stringToDouble(convertData);
	}
	
	double DataConvertor::convertToDouble(string& target, int length)
	{
		length = length < target.size() ? length : target.size();	/* 确保长度不越界 */
		string convertData = target.substr(0, length);/* 获取待转换子串 */
		return stringToDouble(convertData);	/* 转换子串为Int型并返回 */
	}
	/* *************************Function*******************************
	函数名称：convertToString
	函数功能：将从当前转换位置开始，寻找第一个和第二个分割字符, 将中间字符串返回
	输入参数：分割字符的ASCII码
	作者：GDZ
	日期：2021年7月14日 */
	string DataConvertor::convertToString(int length)
	{
		__check(length);
		string ans = (*m_dataPtr).substr(m_transPtr, length);
		m_transPtr += length;
		return ans;
	}
	/* *************************Function*******************************
	函数名称：convertToString
	函数功能：将从当前转换位置开始，寻找第一个和第二个分割字符,将中间字符串返回
	输入参数：分割字符的ASCII码
	作者：GDZ
	日期：2021年7月14日 */
	string DataConvertor::convertToString(int a, int divideSymbsASCII)
	{
		string ans();
		int start = 0;//开始位置
		int length = 0;//子串长度
		while (m_transPtr < m_dataPtr->size()) {
			if ((*m_dataPtr)[m_transPtr] == divideSymbsASCII) {
				start = ++m_transPtr;//开始位置为第一个分割字符位置+1
				break;
			}
			++m_transPtr;
		}
		while (m_transPtr < m_dataPtr->size()) {
			if ((*m_dataPtr)[m_transPtr] == divideSymbsASCII) {
				++m_transPtr;
				break;
			}
			++m_transPtr;
			++length;
		}
		return m_dataPtr->substr(start, length);
	}
	/* *************************Function*************************
	函数名称：__check
	函数功能：检查即将转换的长度是否会超出数据行范围，如果是的那就其缩短
	输入参数：即将转换的长度length
	作者：GDZ
	日期：2021年7月15日 */
	void DataConvertor::__check(int& length)
	{
		length = (m_transPtr + length) <= m_dataPtr->size() ? length : (m_dataPtr->size() - m_transPtr);
	}
	/* *************************Function*************************
	函数名称：stringToInt
	函数功能：将从字符串转换为Int
	输入参数：待转字符串
	作者：GDZ
	日期：2021年7月14日 */
	int DataConvertor::stringToInt(string& target)
	{
		int ans = 0;
		stringstream stream;
		stream << target;
		stream >> ans;
		return ans;
	}

	/* *************************Function*************************
	函数名称：stringToDouble
	函数功能：将从字符串转换为Double
	输入参数：待转字符串
	作者：GDZ
	日期：2021年7月14日 */
	double DataConvertor::stringToDouble(string& target)
	{
		double ans;
		int firstNinusPos = target.find('-');//未找到返回
		int lastNinusPos = target.rfind('-');//未找到返回-1  从后面寻找负号
		int plusPos = target.rfind('+');	//未找到返回-1
		if (firstNinusPos != -1 && firstNinusPos == lastNinusPos) {//数据中只找到一个负号
			if (lastNinusPos != 0 && target.at(lastNinusPos - 1) != ' ') {
				//唯一的负号不在首位          且负号前不是空格     那么说明这个负号是指数的负号			
				double a;	int n;		//数据形式为a×-10^n  a为正数
				stringstream stream_a, stream_n;
				stream_a << target.substr(0, lastNinusPos);
				stream_a >> a;
				stream_n << target.substr(lastNinusPos + 1, target.size() - lastNinusPos - 1);
				stream_n >> n;
				ans = a / (pow(10, n));
			}
			else if (plusPos != -1) {//负号是在第一位 仍有可能后面还有正号
				double a;	int n;		//数据形式为a×10^n  a为负数
				stringstream stream_a, stream_n;
				stream_a << target.substr(0, plusPos);
				stream_a >> a;
				stream_n << target.substr(plusPos + 1, target.size() - plusPos - 1);
				stream_n >> n;
				ans = a / (pow(10, n));
			}
			else if (plusPos == -1) {//负号在第一位
				stringstream stream;
				stream << target;
				stream >> ans;
			}
		}
		else if (firstNinusPos != -1 && firstNinusPos != lastNinusPos) {//找到了位置不同的两个负号
			double a;	int n;		//数据形式为a×-10^n  a为复数
			stringstream stream_a, stream_n;
			stream_a << target.substr(0, lastNinusPos);
			stream_a >> a;
			stream_n << target.substr(lastNinusPos + 1, target.size() - lastNinusPos - 1);
			stream_n >> n;
			ans = a / (pow(10, n));
		}
		else if (plusPos != -1) {//只有一个正号
			double a;	int n;		//数据形式为a×10^n
			stringstream stream_a, stream_n;
			stream_a << target.substr(0, plusPos);
			stream_a >> a;
			stream_n << target.substr(plusPos, target.size() - plusPos - 1);
			stream_n >> n;
			ans = a * (pow(10, n));
		}
		else if (plusPos == -1) {//什么符号都没有
			stringstream stream;
			stream << target;
			stream >> ans;
		}
		return ans;
	}
	;
}