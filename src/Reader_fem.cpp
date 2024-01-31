#include "include/Reader_fem.h"
#include "include/DataConvertor.h"
#include "include/data_management.h"
using namespace std;
namespace CAE 
{
	void Reader_fem::readInputFile(data_management& data_cae)
	{
		this->m_dataConvertor = new DataConvertor();

		this->readData(path_);/* 读入所有数据行 */
		this->classifyData(data_cae);/* 识别数据 */

		data_cae.ne_ = data_cae.node_topos_.size();
		data_cae.nd_ = data_cae.coords_.size();

		string temp_ele_type_str;
		for (std::map<string, int>::iterator it = data_cae.ELE_TYPES.begin(); it != data_cae.ELE_TYPES.end(); it++)
        {
            temp_ele_type_str = it->first;
			data_cae.add_ele(temp_ele_type_str);
        }
        

	}

	/*  函数名称：readData
	函数功能：读文件到fileData，每行保存一个string，忽略"$$""$\node""\node"开头的行
	输入参数：文件路径
	作者：GDZ
	日期：2021年4月14日 */
	void Reader_fem::readData(string file)
	{
		FILE* _targetFile;
		if ((_targetFile = fopen(file.data(), "rt")) == NULL) {
			return;//TODO 异常:打开文件失败
		}
		while (!feof(_targetFile)) {
			char c_line[256];
			fgets(c_line, 255, _targetFile);	//获取行
			if(thisLineIsUseful(c_line)) {//识别是否有效
				string line(c_line);
				line.pop_back();//删除"\node"
				if (line[0] == '+') {			 //此行是否以"+"开始
					this->fileData[fileData.size()-1].append(line.substr(8,line.size()-8));			 //此行添加至上一行
				}
				else {
					this->fileData.push_back(line);//插入行
				}
			}
		}
		fclose(_targetFile);
	}
	bool Reader_fem::thisLineIsUseful(char* line)
	{
		if (strncmp(line, "$$", 2) == 0) { return false; }
		if (strncmp(line, "$\n", 2) == 0) { return false; }
		if (strncmp(line, "\n", 1) == 0) { return false; }

		return true; /* 判断行是否为有效数据 */  
	}
	
	/*  函数名称：classifyData
	函数功能：识别行首关键字，根据关键字调用函数获取数据
	作者：GDZ
	日期：2021年4月16日 */
	void Reader_fem::classifyData(data_management& data_cae)
	{
		for (auto line = this->fileData.begin(); line != this->fileData.end(); line++) {//历遍每一行数据
			if (line->compare(0, sizeof("DESOBJ") - 1, "DESOBJ") == 0) {
				//get_DESOBJ_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("SUBCASE") - 1, "SUBCASE") == 0) {
				//get_SUBCASE_from(line, data_cae);
			}
			if (line->compare(0, sizeof("DTPL") - 1, "DTPL") == 0) {
				//get_DTPL_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("DRESP1") - 1, "DRESP1") == 0) {
				//get_DRESP1_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("DCONSTR") - 1, "DCONSTR") == 0) {
				//get_DCONSTR_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("DCONADD") - 1, "DCONADD") == 0) {
				//get_DCONADD_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("GRID") - 1, "GRID") == 0) {
				get_GRID_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("RBE2") - 1, "RBE2") == 0) {
				//get_RBE2_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("CHEXA") - 1, "CHEXA") == 0) {
				//get_CHEXA_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("CTETRA") - 1, "CTETRA") == 0) {
				get_CTETRA_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("CQUAD4") - 1, "CQUAD4") == 0) {
				get_CQUAD4_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("CTRIA3") - 1, "CTRIA3") == 0) {
				get_CTRIA3_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("$HMNAME COMP") - 1, "$HMNAME COMP") == 0) {
				//get_COMPONENT_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("PSOLID") - 1, "PSOLID") == 0) {
				//get_PSOLID_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("MAT1") - 1, "MAT1") == 0) {
				//get_MAT1_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("SPC    ") - 1, "SPC    ") == 0) {
				get_SPC_from((*line), data_cae);
			}
			if (line->compare(0, sizeof("FORCE") - 1, "FORCE") == 0) {
				get_FORCE_from((*line), data_cae);
			}
	
		}
	}
	//
	///*  函数名称：get_DESOBJ_from
	//函数功能：数据行分割，获取目标函数数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年7月15日 */
	//void Reader_fem::get_DESOBJ_from(string& line, data_management& data_cae)
	//{
	//
	//	fem_OptiObj* _optiObj = new fem_OptiObj();
	//						 m_dataConvertor.inputDataLine(line);	//将此行交给转换器
	//						 m_dataConvertor.ignore(7);				//忽略行首关键字
	//	_optiObj->optiType = m_dataConvertor.convertToString(string::size_type(3));//直接给3会调用ASCII码的重载
	//						 m_dataConvertor.ignore(2);				//忽略)=这两个字符
	//	_optiObj->respID =	 m_dataConvertor.convertToInt(10);		//响应ID
	//	data_cae.insert(_optiObj);
	//}
	//
	///*  函数名称：get_SUBCASE_from
	//函数功能：数据行分割，获取设计变量数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年7月15日 */
	//void Reader_fem::get_SUBCASE_from(vector<string>::iterator& line, data_management& data_cae)
	//{
	//	fem_Subcase* _subCase = new fem_Subcase();
	//	m_dataConvertor.inputDataLine((*line));	//将此行交给转换器
	//	m_dataConvertor.ignore(8);				//忽略行首关键字
	//	_subCase->ID = m_dataConvertor.convertToInt(8);
	//	line += 2;
	//	if (line->compare(9, sizeof("STATICS") - 1, "STATICS") == 0) {
	//		_subCase->type = 1;
	//	}
	//	else if (line->compare(9, sizeof("MFREQ") - 1, "MFREQ") == 0) {
	//		_subCase->type = 2;
	//	}
	//	line++;
	//	m_dataConvertor.inputDataLine((*line));	//将此行交给转换器
	//	m_dataConvertor.ignore(8);				//忽略行首关键字
	//	_subCase->spcID = m_dataConvertor.convertToInt(8);
	//	line++;
	//	m_dataConvertor.inputDataLine((*line));	//将此行交给转换器
	//	m_dataConvertor.ignore(9);				//忽略行首关键字
	//	_subCase->loadID = m_dataConvertor.convertToInt(8);
	//	data_cae.insert(_subCase);
	//}
	//
	///*  函数名称：get_DTPL_from
	//函数功能：获取设计变量数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年7月15日 */
	//void Reader_fem::get_DTPL_from(string& line, data_management& data_cae)
	//{
	//	fem_OptiVar*  _optiVar = new fem_OptiVar();
	//					     m_dataConvertor.inputDataLine(line);	//将此行交给转换器
	//					     m_dataConvertor.ignore(8);				//忽略行首关键字
	//	_optiVar->ID =	     m_dataConvertor.convertToInt(8);
	//	_optiVar->propType = m_dataConvertor.convertToString(string::size_type(8));
	//	_optiVar->porpID =   m_dataConvertor.convertToInt(8);		//属性ID
	//	data_cae.insert(_optiVar);
	//}
	//
	///*  函数名称：get_DRESP1_from
	// 函数功能：获取响应数据
	// 输入参数：数据行
	// 作者：GDZ
	// 日期：2021年7月15日 */
	//void Reader_fem::get_DRESP1_from(string& line, data_management& data_cae)
	//{
	//	fem_OptiResp* _optiResp = new fem_OptiResp();
	//						  m_dataConvertor.inputDataLine(line);	//将此行交给转换器
	//						  m_dataConvertor.ignore(8);			//忽略行首关键字
	//	_optiResp->ID =		  m_dataConvertor.convertToInt(8);
	//	_optiResp->name =	  m_dataConvertor.convertToString(string::size_type(8));
	//	_optiResp->respType = m_dataConvertor.convertToString(string::size_type(8));
	//	data_cae.insert(_optiResp);
	//}
	//
	///*  函数名称：get_DCONSTR_from
	// 函数功能：获取优化约束条件数据
	// 输入参数：数据行
	// 作者：GDZ
	// 日期：2021年7月15日 */
	//void Reader_fem::get_DCONSTR_from(string& line, data_management& data_cae)
	//{
	//	fem_OptiCons* _optiCons = new fem_OptiCons();
	//							m_dataConvertor.inputDataLine(line);//将此行交给转换器
	//							m_dataConvertor.ignore(8);		   //忽略行首关键字
	//	_optiCons->ID =		    m_dataConvertor.convertToInt(8);		 //约束
	//	_optiCons->respID =	    m_dataConvertor.convertToInt(8);		 //响应ID
	//	_optiCons->lowerBound = m_dataConvertor.convertToDouble(8);  //约束下界TODO
	//	_optiCons->upperBound = m_dataConvertor.convertToDouble(8);  //约束上界TODO
	//	data_cae.insert(_optiCons);
	//}
	//
	///*  函数名称：get_DCONADD_from
	// 函数功能：获取约束施加到载荷步数据
	// 输入参数：数据行
	// 作者：GDZ
	// 日期：2021年7月15日 */
	//void Reader_fem::get_DCONADD_from(string& line, data_management& data_cae)
	//{
	//	fem_ConsAdd* _consAdd = new fem_ConsAdd();
	//					   m_dataConvertor.inputDataLine(line);//将此行交给转换器
	//					   m_dataConvertor.ignore(8);//忽略行首关键字
	//	_consAdd->ID =	   m_dataConvertor.convertToInt(8);//ID
	//	_consAdd->consID = m_dataConvertor.convertToInt(8);//约束ID
	//	data_cae.insert(_consAdd);
	//}
	//
	/*  函数名称：get_GRID_from
	 函数功能：获取节点数据
	 输入参数：数据行
	 作者：GDZ
	 日期：2021年4月16日 */
	void Reader_fem::get_GRID_from(string& line, data_management& data_cae)
	{
					m_dataConvertor->inputDataLine(line);	//将此行交给转换器
					m_dataConvertor->ignore(8);				//忽略行首关键字
		int ID =	m_dataConvertor->convertToInt(8);		//节点ID
					m_dataConvertor->ignore(8);				//8个空格
		double x =	m_dataConvertor->convertToDouble(8);	//X
		double y =	m_dataConvertor->convertToDouble(8);	//y
		double z =	m_dataConvertor->convertToDouble(8);	//z

		
		data_cae.coords_.push_back({x,y,z});
	}
	
	///*  函数名称：get_RBE2_from
	// 函数功能：获取RBE2单元数据
	// 输入参数：数据行
	// 作者：GDZ
	// 日期：2021年7月15日 */
	//void Reader_fem::get_RBE2_from(string& line, data_management& data_cae)
	//{
	//	fem_Rbe2* _rbe2 = new fem_Rbe2();
	//						  m_dataConvertor.inputDataLine(line);	//数据行交给转换器
	//						  m_dataConvertor.ignore(8);			//忽略行首关键字
	//	_rbe2->ID =			  m_dataConvertor.convertToInt(8);		//RBE2自己的单元编号
	//	_rbe2->MasterNodeID = m_dataConvertor.convertToInt(8);		//主节点
	//	_rbe2->dof =		  m_dataConvertor.convertToInt(8);		//刚性连接自由度
	//	_rbe2->SlaveNodeID =  m_dataConvertor.convertToInt(8);		//从节点
	//	data_cae.insert(_rbe2);
	//}
	//
	///*  函数名称：get_CHEXA_from
	//函数功能：获取单元数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年4月18日 */
	//void Reader_fem::get_CHEXA_from(string& line, data_management& data_cae)
	//{
	//	fem_Chexa* _chexa = new fem_Chexa();
	//						m_dataConvertor.inputDataLine(line);	//数据行交给转换器
	//						m_dataConvertor.ignore(8);				//忽略行首关键字
	//	_chexa->ID = 		m_dataConvertor.convertToInt(8);		//单元ID
	//	_chexa->propID =	m_dataConvertor.convertToInt(8);		//属性ID
	//	_chexa->nodeID[0] =	m_dataConvertor.convertToInt(8);		//node1_ID
	//	_chexa->nodeID[1] =	m_dataConvertor.convertToInt(8);		//node2_ID
	//	_chexa->nodeID[2] =	m_dataConvertor.convertToInt(8);		//node3_ID
	//	_chexa->nodeID[3] =	m_dataConvertor.convertToInt(8);		//node4_ID
	//	_chexa->nodeID[4] =	m_dataConvertor.convertToInt(8);		//node5_ID
	//	_chexa->nodeID[5] =	m_dataConvertor.convertToInt(8);		//node6_ID
	//	_chexa->nodeID[6] =	m_dataConvertor.convertToInt(8);		//node7_ID
	//	_chexa->nodeID[7] =	m_dataConvertor.convertToInt(8);		//node8_ID
	//	data_cae.insert(_chexa);
	//}
	/*  函数名称：get_CTETRA_from
	函数功能：获取单元数据
	输入参数：数据行
	作者：GDZ
	日期：2021年4月18日 */
	void Reader_fem::get_CTETRA_from(string& line, data_management& data_cae)
	{
		
					  m_dataConvertor->inputDataLine(line);	//数据行交给转换器
					  m_dataConvertor->ignore(8);			//忽略行首关键字
		int ID =	  m_dataConvertor->convertToInt(8);		//单元ID
		int propID  = m_dataConvertor->convertToInt(8);		//属性ID
		int nodeID0 = m_dataConvertor->convertToInt(8);		//node1_ID
		int nodeID1 = m_dataConvertor->convertToInt(8);		//node2_ID
		int nodeID2 = m_dataConvertor->convertToInt(8);		//node3_ID
		int nodeID3 = m_dataConvertor->convertToInt(8);		//node4_ID

	
		data_cae.node_topos_.push_back({ nodeID0 ,nodeID1 ,nodeID2 ,nodeID3 });


		data_cae.ele_list_idx_.push_back(0);
	}
	/*  函数名称：get_CQUAD4_from
	函数功能：获取单元数据
	输入参数：数据行
	作者：GDZ
	日期：2024年1月4日 */
	void Reader_fem::get_CQUAD4_from(string& line, data_management& data_cae)
	{

					  m_dataConvertor->inputDataLine(line);		//数据行交给转换器
					  m_dataConvertor->ignore(8);				//忽略行首关键字
		int ID =	  m_dataConvertor->convertToInt(8);			//单元ID
		int propID =  m_dataConvertor->convertToInt(8);			//属性ID
		int nodeID0 = m_dataConvertor->convertToInt(8);			//node1_ID
		int nodeID1 = m_dataConvertor->convertToInt(8);			//node2_ID
		int nodeID2 = m_dataConvertor->convertToInt(8);			//node3_ID
		int nodeID3 = m_dataConvertor->convertToInt(8);			//node4_ID
		
		data_cae.node_topos_.push_back({ nodeID0 ,nodeID1 ,nodeID2 ,nodeID3 });
	}
	
	/*  函数名称：get_CTRIA3_from
	函数功能：获取单元数据
	输入参数：数据行
	作者：GDZ
	日期：2024年1月4日 */
	void Reader_fem::get_CTRIA3_from(string& line, data_management& data_cae)
	{

		m_dataConvertor->inputDataLine(line);					//数据行交给转换器
		m_dataConvertor->ignore(8);								//忽略行首关键字
		int ID = m_dataConvertor->convertToInt(8);				//单元ID
		int propID = m_dataConvertor->convertToInt(8);			//属性ID
		int nodeID0 = m_dataConvertor->convertToInt(8);			//node1_ID
		int nodeID1 = m_dataConvertor->convertToInt(8);			//node2_ID
		int nodeID2 = m_dataConvertor->convertToInt(8);			//node3_ID

		data_cae.node_topos_.push_back({ nodeID0 ,nodeID1 ,nodeID2});

		
	}
	///* 函数名称：get_Component_from
	//函数功能：获取Component数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年6月25日 */
	//void Reader_fem::get_COMPONENT_from(string& line, data_management& data_cae)
	//{
	//	fem_Component* _component = new fem_Component();
	//						   m_dataConvertor.inputDataLine(line);	//数据行交给转换器
	//						   m_dataConvertor.ignore(24);			//忽略行首关键字	
	//	_component->ID =	   m_dataConvertor.convertToInt(8);		//单元组ID
	//	_component->name =     m_dataConvertor.convertToString(34);	//组件名称 组件名称在文件中是用""分割，因此输入为"的ASCII码34
	//						   m_dataConvertor.ignore(1);			//忽略fem文件名称后面的一个多余空格
	//	_component->propID =   m_dataConvertor.convertToInt(8);		//属性ID
	//	_component->propName = m_dataConvertor.convertToString(34);	//属性名称
	//						   m_dataConvertor.ignore(1);			//忽略fem文件名称后面的一个多余空格
	//						   m_dataConvertor.convertToInt(2);		//不知名ID
	//	data_cae.insert(_component);
	//}
	//
	///* 函数名称：get_PSOLID_from
	//函数功能：获取实体属性数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年7月16日 */
	//void Reader_fem::get_PSOLID_from(string& line, data_management& data_cae)
	//{
	//	fem_Prop* _psolid = new fem_Prop();	
	//					 m_dataConvertor.inputDataLine(line);	//数据行交给转换器
	//					 m_dataConvertor.ignore(8);			//忽略行首关键字	
	//	_psolid->ID =	 m_dataConvertor.convertToInt(8);		//单元组ID
	//	_psolid->matID = m_dataConvertor.convertToInt(8);		//属性ID
	//	data_cae.insert(_psolid);
	//}
	//
	///*  函数名称：get_MAT1_from
	//函数功能：获取材料数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年7月16日 */
	//void Reader_fem::get_MAT1_from(string& line, data_management& data_cae)
	//{
	//	fem_Mat1* _mat1 = new fem_Mat1();
	//			m_dataConvertor.inputDataLine(line);	//数据行交给转换器
	//			m_dataConvertor.ignore(8);			//忽略行首关键字	
	//	_mat1->ID =  m_dataConvertor.convertToInt(8);		//单元组ID
	//	_mat1->E =   m_dataConvertor.convertToDouble(8);
	//	_mat1->G =   m_dataConvertor.convertToDouble(8);
	//	_mat1->NU =  m_dataConvertor.convertToDouble(8);
	//	_mat1->RHO = m_dataConvertor.convertToDouble(8);
	//	_mat1->A =	 m_dataConvertor.convertToDouble(8);
	//	_mat1->TREF= m_dataConvertor.convertToDouble(8);
	//	_mat1->GE =  m_dataConvertor.convertToDouble(8);
	//	_mat1->ST =  m_dataConvertor.convertToDouble(8);
	//	_mat1->SC =  m_dataConvertor.convertToDouble(8);
	//	_mat1->SS =  m_dataConvertor.convertToDouble(8);
	//	data_cae.insert(_mat1);
	//}
	//
	///*  函数名称：get_SPC_from
	//函数功能：获取单点约束数据
	//输入参数：数据行
	//作者：GDZ
	//日期：2021年4月26日 */
	void Reader_fem::get_SPC_from(string& line, data_management& data_cae)
	{
		
					   m_dataConvertor->inputDataLine(line);	//数据行交给转换器
					   m_dataConvertor->ignore(8);				//忽略行首关键字
		int colID =    m_dataConvertor->convertToInt(8);		//所属Load Collector ID
		int nodeID =   m_dataConvertor->convertToInt(8);		//作用节点ID
		int _dof =	   m_dataConvertor->convertToInt(8);		//约束自由度	
		
		data_cae.dis_bc_set_.push_back(nodeID);
	}


	void Reader_fem::transSpcDof(bool spcDof[6], int _dof)
	{
		while (_dof != 0) {
			spcDof[_dof % 10-1] = true;
			_dof = _dof / 10;
		}
	}
	
	/*  函数名称：get_FORCE_from
	函数功能：获取载荷数据
	输入参数：数据行
	作者：GDZ
	日期：2021年4月26日  */
	void Reader_fem::get_FORCE_from(string& line, data_management& data_cae)
	{
	
					     m_dataConvertor->inputDataLine(line);		//数据行交给转换器
					     m_dataConvertor->ignore(8);				//忽略行首关键字
		int colID =      m_dataConvertor->convertToInt(8);			//所属Load Collector ID
		int nodeID =     m_dataConvertor->convertToInt(8);			//作用节点ID
		 			     m_dataConvertor->convertToInt(8);			//坐标系编号	**暂时用不上
		 			     m_dataConvertor->convertToDouble(8);		//比例因子	**暂时用不上
		double xComp =	 m_dataConvertor->convertToDouble(8);		//x方向分量
		double yComp =	 m_dataConvertor->convertToDouble(8);		//y方向分量
		double zComp =	 m_dataConvertor->convertToDouble(8);		//z方向分量

		data_cae.load_set_.push_back(nodeID);

		if (abs(xComp) > 0.00001) {
			data_cae.load_dof_ = 1;
			data_cae.load_value_ = xComp;
		}
		else if (abs(yComp) > 0.00001) {
			data_cae.load_dof_ = 2;
			data_cae.load_value_ = yComp;
		}
		else if (abs(zComp) > 0.00001) {
			data_cae.load_dof_ = 3;
			data_cae.load_value_ = zComp;
		}
	}
}