# cmake最低版本号要求
cmake_minimum_required(VERSION 3.26)

# 设置项目名
set(PROJECT_NAME CAE)

# 添加src文件夹下的所有源文件
aux_source_directory(${CMAKE_SOURCE_DIR}/src SRC_FILE)

# 设置可执行文件储存路径
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

# 编译
add_executable(${PROJECT_NAME} ${SRC_FILE})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/eigen-3.4.0)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/externals/superlu-6.0.1)
target_link_libraries(${PROJECT_NAME} PUBLIC superlu)
target_link_libraries(${PROJECT_NAME} PUBLIC blas)
