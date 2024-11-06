# cmake最低版本号要求
cmake_minimum_required(VERSION 3.23)

# 使VS目录和文件夹目录一致
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

# 添加头文件查找路径
include_directories(./)
include_directories(./Elements)
include_directories(./solver)
include_directories(./include)
include_directories(./Elements/include)
include_directories(./solver/include)

# 添加src文件夹下的所有源文件
file(GLOB_RECURSE SRC_FILE CONFIGURE_DEPENDS *.cpp *.h)
source_group_by_dir(SRC_FILE)# 调用函数使VS目录和文件夹目录一致

# 设置可执行文件储存路径
SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin/")
if(MSVC)
    #在Visual Studio生成器（即VS工程）下，会在EXECUTABLE_OUTPUT_PATH、EXECUTABLE_OUTPUT_PATH后面自动加一个${CMAKE_BUILD_TYPE}
    SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/) #可执行文件输出目录
    set(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/externals)    #库文件输出目录
else()
    #为了和VS相同，在其他环境下，EXECUTABLE_OUTPUT_PATH设置为${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}
    set(CMAKE_INSTALL_RPATH "$ORIGIN/externals/mkl:$ORIGIN/externals")#库文件相对路径，仅linux有这个功能
    SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
    set(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/externals)
endif()

# 编译可执行程序
add_executable(${PROJECT_NAME} ${SRC_FILE})

# 链接 MKL
if(MKL_USE)
    target_link_libraries(${PROJECT_NAME} PUBLIC MKL::MKL)
endif()

# 链接 Eigen
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/eigen-3.4.0)

# 链接 SuperLU
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/superlu-6.0.1)
target_link_libraries(${PROJECT_NAME} PUBLIC superlu)
target_link_libraries(${PROJECT_NAME} PUBLIC blas)



# # cmake最低版本号要求
# cmake_minimum_required(VERSION 3.23)

# # 使VS目录和文件夹目录一致
# macro(source_group_by_dir source_files)
#     if(MSVC)
#         set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR})
#         foreach(sgbd_file ${${source_files}})
#             string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
#             string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
#             string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
#             string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
#             if(sgbd_nogroup)
#                 set(sgbd_group_name "\\")
#             endif(sgbd_nogroup)
#             source_group(${sgbd_group_name} FILES ${sgbd_file})
#         endforeach(sgbd_file)
#     endif(MSVC)
# endmacro(source_group_by_dir)

# # 添加头文件查找路径
# include_directories(./)
# include_directories(./Elements)
# include_directories(./solver)
# include_directories(./include)
# include_directories(./Elements/include)
# include_directories(./solver/include)

# # 添加src文件夹下的所有源文件
# file(GLOB_RECURSE SRC_FILE CONFIGURE_DEPENDS *.cpp *.h)
# source_group_by_dir(SRC_FILE)# 使VS目录和文件夹目录一致

# # 设置可执行文件储存路径
# if(MSVC)
#     #在Visual Studio生成器（即VS工程）下，会在EXECUTABLE_OUTPUT_PATH、EXECUTABLE_OUTPUT_PATH后面自动加一个${CMAKE_BUILD_TYPE}
#     SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/) #可执行文件输出目录
#     set(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/externals)    #库文件输出目录
# else()
#     #为了和VS相同，在其他环境下，EXECUTABLE_OUTPUT_PATH设置为${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}
#     set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
#     set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#     set(CMAKE_INSTALL_RPATH "$ORIGIN/externals")#库文件相对路径-全局
#     SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
#     set(LIBRARY_OUTPUT_PATH ${EXECUTABLE_OUTPUT_PATH}/externals)
# endif()

# # 编译可执行程序
# add_executable(${PROJECT_NAME} ${SRC_FILE})

# # 链接 MKL
# if(MSVC)
#     SET(MKL_LIBDIR "${PROJECT_SOURCE_DIR}/externals/mkl/lib")                        
#     file(GLOB MKL_LIBS "${MKL_LIBDIR}/*.lib")
#     SET(MKL_INC "${PROJECT_SOURCE_DIR}/externals/mkl/include")        
#     target_include_directories(${PROJECT_NAME} PUBLIC ${MKL_INC})
#     foreach(MKL_LIB ${MKL_LIBS})
#         target_link_libraries(${PROJECT_NAME} PUBLIC ${MKL_LIB})
#     endforeach()
# else()
#     set_target_properties(MKL::MKL PROPERTIES INSTALL_RPATH "$ORIGIN/externals/mkl")#运行时mkl库文件查找路径
#     target_link_libraries(${PROJECT_NAME}  PUBLIC MKL::MKL)
#     # 将MKL的库拷贝到：可执行程序路径/externals/mkl
#     foreach(mkl_ IN LISTS MKL::MKL)
#         get_filename_component(mkl_lib_name "${mkl_}" NAME)
#         file(COPY "${mkl_lib_name}" DESTINATION "${LIBRARY_OUTPUT_PATH}/mkl")
#     endforeach()
# endif(MSVC)

# # 链接 Eigen
# target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/eigen-3.4.0)

# # 链接 SuperLU
# target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/superlu-6.0.1)
# target_link_libraries(${PROJECT_NAME} PUBLIC superlu)
# target_link_libraries(${PROJECT_NAME} PUBLIC blas)