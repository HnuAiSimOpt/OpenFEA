# cmake最低版本号要求
cmake_minimum_required (VERSION 3.23)

# 设置项目名
set(PROJECT_NAME OpenFEA)
project(${PROJECT_NAME})

# 设置编译器参数
set(CMAKE_SKIP_RPATH false)
set(CMAKE_SKIP_BUILD_RPATH false)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set (CMAKE_CXX_STANDARD_REQUIRED ON)  
set (CMAKE_CXX_STAND 17)

# 添加superlu子文件
set(enable_internal_blaslib ON)# 设置superLU使用内部blas库，避免与mkl的冲突
add_subdirectory(externals/superlu-6.0.1)

# 查找mkl
set(MKL_USE TRUE)
if(MKL_USE)
    list(APPEND CMAKE_PREFIX_PATH $ENV{MKLROOT}/lib/cmake/mkl)
    set (MKL_INTERFACE lp64)
    set (MKL_LINK static) # 使用静态链接，生成的程序不用拷贝动态库（如果希望使用动态链接，请查阅CmakeLists参数手册、打包发送注意事项）
    # set (MKL_THREADING tbb_thread)
    find_package(MKL CONFIG REQUIRED)
endif()


add_subdirectory(src)



# # cmake最低版本号要求
# cmake_minimum_required (VERSION 3.23)

# SET(CMAKE_C_COMPILER_WORKS TRUE)
 
# SET(CMAKE_CXX_COMPILER_WORKS TRUE)
# # 设置项目名
# set(PROJECT_NAME OpenFEA)
# project(${PROJECT_NAME})
# SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin/unpack")



# SET(CMAKE_SKIP_RPATH false)
# SET(CMAKE_SKIP_BUILD_RPATH false)
# SET(CMAKE_BUILD_RPATH_USE_ORIGIN True)
# # 设置编译器参数
# set (CMAKE_CXX_STANDARD_REQUIRED ON)  
# set (CMAKE_CXX_STAND 17)


# # 添加子文件
# set(enable_internal_blaslib ON)# 设置superLU使用内部blas库，避免与mkl的冲突
# add_subdirectory(externals/superlu-6.0.1)

# if(MSVC)
#     add_subdirectory(externals/mkl)
# else()
#     list(APPEND CMAKE_PREFIX_PATH $ENV{MKLROOT}/lib/cmake/mkl)
#     set (MKL_INTERFACE lp64)
#     #set (MKL_LINK Static)
#     find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
# endif(MSVC)

# add_subdirectory(src)
