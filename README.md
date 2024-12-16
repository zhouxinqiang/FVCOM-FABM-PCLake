# FVCOM-FABM-PCLake
To build FVCOM-FABM-PCLake coupled model, first compile FABM and PCLake and then link FVCOM to it.

## Complite FABM and PCLake

1.Download FABM and PCLake  
2.Extract the FABM source code (the following example code extracts the code to $HOME/Code/fabm/src)  
3.Extract the PCLake source code (the following example code extracts the code to $HOME/Code/PCLake/src)  
4.Ensure that CMake 2.8.8 or later is installed.  
5.Compile FABM and PCLake  

```
  cd $HOME/Code/fabm/src
  mkdir build
  cd build
  cmake $HOME/Code/fabm/src -DFABM_HOST=fvcom -DFABM_PCLake_BASE=$HOME/Code/PCLake -DCMAKE_Fortran_COMPILER=$(which mpif90)
  make install
```

## Complite FVCOM 

Compile FVCOM with FABM support, please edit make.inc as follows to set FLAG_25 and compile FVCOM normally. If FABM has been installed to a custom location, please adjust the BIOLIB and BIOINCS paths as needed.

```
            # Online configuration  
            FLAG_25 = -DFABM  
            BIOLIB       = -L$(HOME)/local/fabm/fvcom/lib -lfabm  
            BIOINCS      = -I$(HOME)/local/fabm/fvcom/include  
```

# 1. FVCOM编译 

编译主要代码及其库：  

打开一个终端。  
设置MPI环境。在Fedora / Red Hat / CentOS上， yum install mpich 提供了一个mpich环境，可以将其加载 module load mpi/mpich-x86_64。  
更改下载FVCOM的目录并解压缩代码。 
输入 FVCOM_source目录。 
编辑  make.inc 文件以启用/禁用其他功能，例如润湿/干燥以满足您的要求。  
将     TOPDIR 变量更改为您当前所在的路径（即 pwdLinux命令的输出）。  
输入   libs 子目录。  
输入  make 并等待编译完成。  
切换回父  FVCOM_source 目录。  
输入  make 并等待编译完成。  
将  fvcom 二进制文件复制到模型输入目录，然后切换到该目录。  
使用以下命令启动模型 mpirun：  
 
    mpirun -n $num_proc fvcom --casename test --dbg=0 --logfile=fvcom.log
$num_proc 处理器的数量 在哪里。

该 make.inc文件包含可能感兴趣的一系列环境变量，特别是编译器选项。参见第463-469行。对于gfortran，以下环境变量可make.inc用于Fedora。

```
#--------------------------------------------------------------------------
#  MPIF90/GFORTRAN Compiler Definitions (PML)
#--------------------------------------------------------------------------
         CPP      = cpp
         COMPILER = -DGFORTRAN
         FC       = mpif90
         DEBFLGS  =
         OPT      = -O3 -L/usr/lib64/mpich-x86_64/lib -I/usr/include/mpich-x86_64 -I/usr/lib64/gfortran/modules/
         CLIB     =
         CC       = mpicc
```

## 2. 非静液压和半隐式FVCOM
Non-hydrostatic and semi-implicit FVCOM

 
### 2.1 依存关系
对于非静水，半隐式，数据同化，卡尔曼滤波器和波电流相互作用，FVCOM需要PETSc和HYPRE。FVCOM是根据PETSc库的2.3.3版本编写的; PETSc版本3.x不适用于FVCOM。反过来，PETSc依赖于HYPRE，并且我们已经成功使用了HYPRE 2.0.0版。

### 2.2 编译
在您的中make.inc，在TOPDIR声明之后添加一个新部分，其中包括两个新变量：

```
#--------------------------------------------------------------------------
#  PETSC library locations (for non-hydrostatic/semi-implicit/data assimilation)
#--------------------------------------------------------------------------
            PETSC_LIB     =  -L$(PETSC_DIR)/lib/linux-gnu-intel/
            PETSC_FC_INCLUDES =  -I$(PETSC_DIR) -I$(PETSC_DIR)/bmake/$(PETSC_ARCH) -I$(PETSC_DIR)/include
```

确保您的PETSc安装将PETSC_DIR和PETSC_ARCH环境变量设置为PETSc安装的根目录（如果您自己编译了PETSc，PETSC_DIR则是使用--prefix in 指定的路径configure.py）。如果PETSC_ARCH未定义，则应该可以通过查找$PETSC_DIR/lib/; 来识别有效值。有效值将是该目录中目录的名称。

然后，对于非流体力学，请编辑您的内容make.inc以取消注释该FLAG_30 = -DNH行及其include $(PETSC_DIR)/bmake/common/variables下方以及该FLAG_9 = -DSEMI_IMPLICIT行。对于其他选项（非静液压，半隐式，数据同化，卡尔曼滤波器和波电流相互作用），请取消注释相关的FLAGs。

最后，将PETSC_LIB和PETSC_FC_INCLUDES变量附加到LIBS和INCLUDES定义的末尾make.inc：

```
            LIBS  = $(LIBDIR) $(CLIB)  $(PARLIB) $(IOLIBS)  $(DTLIBS)\
            $(MPILIB) $(GOTMLIB) $(KFLIB) $(BIOLIB) \
            $(OILIB) $(VISITLIB) $(PROJLIBS) $(PETSC_LIB)
 
            INCS  =     $(INCDIR) $(IOINCS) $(GOTMINCS) $(BIOINCS)\
             $(VISITINCPATH) $(PROJINCS) $(DTINCS) \
             $(PETSC_FC_INCLUDES)
```
要编译FVCOM，请make像往常一样键入。

## 3. FVCOM-GOTM
可在GOTM网站上找到GOTM的安装指南。

这些说明假定您使用的是Intel Fortran编译器（ifort）。

按照上面网站上的说明下载并提取GOTM存档。
打开一个终端窗口，将目录更改为GOTM源代码，然后键入以下内容：
```
    export NETCDFINC=/YOUR/FVCOM/LIBS/DIR/libs/install/include
    export NETCDFLIBNAME=/YOUR/FVCOM/LIBS/DIR/libs/install/lib/libnetcdf.a
    export GOTMDIR=$(pwd)/gotm-4.0.0
    export FORTRAN_COMPILER=IFORT
    cd src
    make
```
用 $TOPDIR 里面的make.inc将FVCOM中/YOUR/FVCOM/LIBS/DIR的值替换。 

这将产生类似于上面网站所述的输出，以及一个名为  gotm_prod_IFORT
通过从GOTM网站下载并运行测试用例来测试您的GOTM构建。
要将GOTM与FVCOM链接，您需要对FVCOM进行以下更改make.inc：
```
    FLAG_11 = -DGOTM
    GOTMLIB       = -L/YOUR/GOTM/DIR/gotm-4.0.0/lib/IFORT/ -lturbulence_prod -lutil_prod -lmeanflow_prod
    GOTMINCS      = -I/YOUR/GOTM/DIR/gotm-4.0.0/modules/IFORT/
```
和以前一样，/YOUR/GOTM/DIR/使用包含GOTM版本的目录进行调整。

在.nml文件中设置以下选项：
```
    BOTTOM_ROUGHNESS_TYPE           = 'gotm'
```
最后，将GOTM输入包含casename_gotmturb.inp在与其他FVCOM输入相同的目录中的文件中。


## 4. FVCOM-FABM：支持生物地球化学模型
要构建FVCOM-FABM，首先要编译FABM，然后将FVCOM链接到它。要下载FVCOM-FABM代码，请先向UMASSD 注册，然后向PML 注册，以访问FVCOM-FABM代码。注册后，FABM-ERSEM从PML GitLab存储库下载分支。

### 4.1 FABM
从https://github.com/fabm-model/fabm下载FABM 。
提取FABM源代码（以下示例代码将代码提取到$HOME/Code/fabm/src）
如果您想将FVCOM-FABM与ERSEM结合使用，请注册访问ERSEM源代码，从PML GitLab下载稳定的代码，然后提取源代码（下面的示例代码将代码提取到$HOME/Code/ersem）
确保已安装CMake 2.8.8或更高版本。
编译FABM：
```
cd $HOME/Code/fabm/src
mkdir build
cd build
cmake $HOME/Code/fabm/src -DFABM_HOST=fvcom -DFABM_ERSEM_BASE=$HOME/Code/ersem -DCMAKE_Fortran_COMPILER=$(which mpif90)
make install
```
-DFABM_ERSEM_BASE=...如果您不使用ERSEM，请 忽略该开关。

要启用调试版本，请将cmake命令更改为以下内容：
```
cmake $HOME/Code/fabm/src -DFABM_HOST=fvcom -DFABM_ERSEM_BASE=$HOME/Code/ersem -DCMAKE_Fortran_COMPILER=$(which mpif90) -DCMAKE_BUILD_TYPE=debug -DCMAKE_Fortran_FLAGS_DEBUG="-g -traceback -check all"

```
在上文中，的值-DCMAKE_Fortran_FLAGS_DEBUG特定于Intel Fortran编译器；如果您使用其他编译器，则将其替换为适合使用该编译器进行调试的标志。

默认情况下，FABM的构建具有双精度。如果-DDOUBLE_PRECISION在编译FVCOM时使用该标志，则这是适当的。如果您打算改用单精度，则需要添加-DFABM_REAL_KIND='SELECTED_REAL_KIND(6)'对的调用cmake。

### 4.2 FVCOM-FABM
要在FABM支持下编译FVCOM，请FLAG_25按如下所示编辑make.inc并正常编译FVCOM。如果已将FABM安装到自定义位置，请根据需要调整BIOLIB和BIOINCS路径。
```
            # Online configuration
            FLAG_25 = -DFABM
            BIOLIB       = -L$(HOME)/local/fabm/fvcom/lib -lfabm
            BIOINCS      = -I$(HOME)/local/fabm/fvcom/include
To enable offline FVCOM-ERSEM runs, change FLAG_25 to -DFABM -DOFFLINE_FABM.
```
