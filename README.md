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
