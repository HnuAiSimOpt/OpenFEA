#pragma once
#include <iostream>
#include <map>

namespace CAE
{
	//使用方法：
	//基类类型* 变量名 = GetFatory<基类类型Factory>::instance()->create_class(key);
	template<typename T>
	class GetFactory
	{
	public:
		static T* instance() {
			if (m_instance == NULL) {
				m_instance = new T();
			}
			return m_instance;
		}
	private:
		GetFactory() {}
		GetFactory(const GetFactory<T>&);
		GetFactory<T>& operator =(const GetFactory<T>&);
		~GetFactory() {}
	private:
		static T* m_instance;
	};

	template<typename T>
	T* GetFactory<T>::m_instance = NULL;


	#define CREAT_FACTORY(ClassBase)											\
	typedef ClassBase * (* create_##ClassBase)(void);							\
	class ClassBase##Factory{													\
		friend class GetFactory<ClassBase##Factory>;							\
	public:																		\
		void register_class(const std::string& className, create_##ClassBase method){	\
			m_classMap[className] = method;									\
		}																\
		ClassBase* create_class(const std::string& className)					\
		{																				\
			auto it = m_classMap.find(className);								\
			if (it != m_classMap.end()) {return it->second();}					\
			return nullptr;														\
		}																		\
		private:																\
			ClassBase##Factory() {}												\
			~ClassBase##Factory() {}											\
		private:																\
			std::map<std::string, create_##ClassBase> m_classMap;								\
	};																			\
	class ClassRegister {												\
		public:																\
			ClassRegister(const std::string& className, create_##ClassBase method)	\
			{																\
				GetFactory<ClassBase##Factory>::instance()->register_class(className, method);\
			}																\
		};

	//在子类的cpp中使用，在h文件中可能会造成重复定义
	#define REGISTER(ClassBase,ClassName,ClassKey)									\
	ClassBase * createObject##ClassName()			\
	{											\
		ClassBase* obj = new ClassName();			\
		return obj;								\
	}											\
	ClassRegister classRegister##ClassBase##ClassName(ClassKey, createObject##ClassName)



}//namespace CAE
