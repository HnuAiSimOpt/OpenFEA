# OpenFEA
## 下载须知
如需完整代码请使用 git clone。  
直接下载代码压缩包，因git尾序LF与CRLF差异可能无法在windows下编译（如确有需求，可自行转换尾序为CRLF）  

## 编译前准备
- cmake
3.23及以上版本
- mkl（可选，建议安装）
安装最新版本Intel oneAPI MKL，并配置相关环境变量
  - *本项目默认使用MKL-Pardiso作为矩阵求解器，如未安装MKL，需修改* {Project}/CMakeLists.txt:24 *，改为*
    > set(MKL_USE FALSE)
  - 将/src/linear_solution.cpp、/src/solver/include/solver_pardiso.h、/src/solver/solver_pardiso.cpp相关内容注释掉

- MSVC（仅Windows）
可直接安装Visual Studio 2019以上版本
- gcc（仅Linux）
11.x及以上版本
  > scl enable devtoolset-11 bash
## 编译
1. 克隆本仓库
   > git clone `https://github.com/HnuAiSimOpt/OpenFEA.git`
2. cmake编译
   **(1) Windows**
      1. **VS Code + Cmake命令**
         - VS code 打开git项目所在目录文件夹
         - 如已安装mkl，设置mkl环境(使用intel提供环境设置脚本，按oneapi实际安装路径修改下行命令)
           在vs code终端运行
           > *{ONEAPI_PATH}*/setvars.bat  

           可使用以下命令检查是否设置成功(成功会返回路径结果)
           > set MKLROOT
         - 新建并进入build文件夹，可在终端输入以下命令
              > mkdir build
              > cd build
         - 使用cmake编译，下面命令为编译一个Release版本的OpenFEA（默认静态链接mkl）
              > cmake ../ -DCMAKE_BUILD_TYPE=Release
              > cmake --build ./ -j *$(nproc)*
      2. **Visual Studio + CmakeGUI**
         - 如已安装mkl，设置mkl环境(使用intel提供环境设置脚本，按oneapi实际安装路径修改下行命令)
           在cmd运行
           > *{ONEAPI_PATH}*/setvars.bat  

           可使用以下命令检查是否设置成功(成功会返回路径结果)
           > set MKLROOT
         - 在cmake-gui设置相关路径及编译器，Configure设置相关编译选项，Generate-Open Project即可

    **(2) Linux**
      - 打开git项目所在目录文件夹
      - 如已安装mkl，设置mkl环境(使用intel提供环境设置脚本，如oneapi安装路径不一致，按实际路径修改下行命令)
          > source /opt/intel/oneapi/setvars.sh   

          可使用以下命令检查是否设置成功(成功会返回路径结果)
          > printenv MKLROOT
      - 新建并进入build文件夹
           > mkdir build
           > cd build
      - 使用cmake编译，下面命令为编译一个Release版本的OpenFEA（默认静态链接mkl）
              > cmake ../ -DCMAKE_BUILD_TYPE=Release
              > cmake --build ./ -j *$(nproc)*
1. 编译完成
cmake编译成功无错后，可在 **{Project}/bin/** 目录下看到编译成功的程序 OpenFEA.exe(Linux无后缀)

## 运行
在终端输入以下命令进行求解
> **OpenFEAPath** -i **InputFilePath** [-o] [**OutputFilePath**] [···]   

**OpenFEAPath**：编译成功的可执行程序路径
**InpFilePath**：求解文件路径（eg: ./input.inp）
**OutputFilePath**：输出结果文件路径（eg：./result.vtk）
*可选命令*
-imp：隐式求解（可不加，默认为隐式计算）
-exp：显式求解
-sfem：光滑有限元
-nl：开启几何非线性计算（默认关闭）