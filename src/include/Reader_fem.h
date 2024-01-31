#pragma once
#include <string>
#include <vector>
#include "data_management.h"
using namespace std;

namespace CAE
{
	class DataConvertor;
	class data_management;
	class Reader_fem
	{

	public:
		Reader_fem(string path) : path_(path) {};
		~Reader_fem() {};

		void readInputFile(data_management& data_cae);
	private:
		string path_;						//�ļ�·��
		vector<string>	fileData;			//�ļ��е�����
		DataConvertor*	m_dataConvertor;	//����ת����

		int C3D4Num = -1;


		void readData(string file);
		bool thisLineIsUseful(char* line);
		void classifyData(data_management& model);
		//void get_DESOBJ_from(string& line, data_management& data_cae);
		//void get_SUBCASE_from(vector<string>::iterator& line, data_management& data_cae);
		//void get_DTPL_from(string& line, data_management& data_cae);
		//void get_DRESP1_from(string& line, data_management& data_cae);
		//void get_DCONSTR_from(string& line, data_management& data_cae);
		//void get_DCONADD_from(string& line, data_management& data_cae);
		void get_GRID_from(string& line, data_management& data_cae);
		//void get_RBE2_from(string& line, data_management& data_cae);
		//void get_CHEXA_from(string& line, data_management& data_cae);
		void get_CTETRA_from(string& line, data_management& data_cae);
		void get_CQUAD4_from(string& line, data_management& data_cae);
		void get_CTRIA3_from(string& line, data_management& data_cae);
		//void get_COMPONENT_from(string& line, data_management& data_cae);
		//void get_PSOLID_from(string& line, data_management& data_cae);
		//void get_MAT1_from(string& line, data_management& data_cae);
		void get_SPC_from(string& line, data_management& data_cae);
		void transSpcDof(bool spcDof[6], int _dof);
		void get_FORCE_from(string& line, data_management& data_cae);
	};
}

