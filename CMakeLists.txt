# cmake最低版本号要求
cmake_minimum_required (VERSION 3.5)

# 设置项目名
set(PROJECT_NAME project_main)
project(${PROJECT_NAME})

# 设置编译器参数
set (CMAKE_CXX_STAND 17)
set (CMAKE_CXX_STANDARD_REQUIRED OFF)  

# 添加子文件
add_subdirectory(src)
add_subdirectory(externals/mkl)
add_subdirectory(externals/superlu-6.0.1)
