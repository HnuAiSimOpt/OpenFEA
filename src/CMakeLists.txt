# cmake最低版本号要求
cmake_minimum_required(VERSION 3.5)

# 设置项目名
set(PROJECT_NAME CAE)

# 使vs目录和文件夹目录一致
macro(source_group_by_dir source_files)
    if(MSVC)
        set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR})
        foreach(sgbd_file ${${source_files}})
            string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
            string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
            string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
            string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
            if(sgbd_nogroup)
                set(sgbd_group_name "\\")
            endif(sgbd_nogroup)
            source_group(${sgbd_group_name} FILES ${sgbd_file})
        endforeach(sgbd_file)
    endif(MSVC)
endmacro(source_group_by_dir)

# 添加src文件夹下的所有源文件
#aux_source_directory(${CMAKE_SOURCE_DIR}/src SRC_FILE)
file(GLOB_RECURSE SRC_FILE CONFIGURE_DEPENDS *.cpp *.h)
source_group_by_dir(SRC_FILE)

# 设置可执行文件储存路径
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

# 编译
add_executable(${PROJECT_NAME} ${SRC_FILE})

# 链接 MKL
SET(MKL_LIBDIR "${PROJECT_SOURCE_DIR}/externals/mkl/lib")                        
file(GLOB MKL_LIBS "${MKL_LIBDIR}/*.lib")
SET(MKL_INC "${PROJECT_SOURCE_DIR}/externals/mkl/include")        
target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INC})
foreach(MKL_LIB ${MKL_LIBS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${MKL_LIB})
endforeach()
# 链接 Eigen
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/eigen-3.4.0)
# 链接 SuperLU
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/superlu-6.0.1)
target_link_libraries(${PROJECT_NAME} PUBLIC superlu)
target_link_libraries(${PROJECT_NAME} PUBLIC blas)