############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   May 10, 2015
#
#  Modified:	    
#		    
#
############################################################################
#
#  The name of the libraries to be created/linked to
#
# SuperLUroot	= E:/CADCAE_project/OpenFEA/build
#
SuperLUroot	= C:/Program Files/project_main
SUPERLULIB   	= $(SuperLUroot)/SRC/libsuperlu.a

TMGLIB       	= libtmglib.a

XSDK_INDEX_SIZE = 
HAVE_METIS      = 

# BLASDEF 	= -DUSE_VENDOR_BLAS
BLASLIB		= C:/Program Files/project_main/CBLAS/libblas.a
LIBS		= $(SUPERLULIB) C:/Program Files/project_main/CBLAS/libblas.a
LIBS  += 

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = D:/vs2019/VC/Tools/MSVC/14.29.30133/bin/Hostx64/x64/lib.exe
ARCHFLAGS    = cr
RANLIB       = :

CC           = D:/vs2019/VC/Tools/MSVC/14.29.30133/bin/Hostx64/x64/cl.exe
CFLAGS 	     = /DWIN32 /D_WINDOWS 
NOOPTS       = -O0
FORTRAN	     = 

LOADER       = $(CC)
LOADOPTS     =

#
# The directory in which Matlab is installed
#
MATLAB	     = /Applications/MATLAB_R2015b.app
