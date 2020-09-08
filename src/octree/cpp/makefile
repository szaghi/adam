MACHINE = galileo
ENABLE_CUDA = 0
TARGET = miniamr
OPENMP = -qopenmp
ifeq ($(ENABLE_CUDA),1)
    #CMD="g++ -O3 -std=c++11"
    #CMD="g++ -g -O0 -C -std=c++11"
    CXX = icpc
    CUDACXXFLAGS = -DUSE_CUDA -DNVIDIA_GPU -DDEBALLOC
    CXXFLAGS = -O3 -std=c++11 -g $(OPENMP) $(CUDACXXFLAGS)
    CC = icc
    NVCCFLAGS = -arch=sm_35
    NVCCFLAGS = -arch=sm_70
    #CUDALDFLAGS = -L/cineca/prod/compilers/cuda/9.0/none/lib64 -lcudart 
    #CUDAINCLUDE = -I/cineca/prod/compilers/cuda/9.0/none/include
    CUDALDFLAGS = -L/cineca/prod/compilers/cuda/8.0.61/none/lib64 -lcudart 
    CUDAINCLUDE = -I/cineca/prod/compilers/cuda/8.0.61/none/include
    CUDA_PATH = /cineca/prod/compilers/cuda/8.0.61/none/
    CUDALDFLAGS = -L/cineca/prod/compilers/cuda/9.0/none/lib64 -lcudart 
    CUDAINCLUDE = -I/cineca/prod/compilers/cuda/9.0/none/include
    CUDA_PATH = /cineca/prod/compilers/cuda/9.0/none/
    LD = icpc -std=c++11 $(CUDALDFLAGS)
    OPENMP = 
    ifeq ($(MACHINE),davide)
        CXX = clang++
        CC = clang
        CUDALDFLAGS = -L/cineca/prod/opt/compilers/cuda/9.2.88/none/lib64 -lcudart 
        CUDAINCLUDE = -I/cineca/prod/opt/compilers/cuda/9.2.88/none/include
        CUDA_PATH = /cineca/prod/opt/compilers/cuda/9.2.88/none/
        LD = clang++ -std=c++11 $(CUDALDFLAGS)
        NVCCFLAGS = -arch=sm_60
    endif
else ifeq ($(ENABLE_CUDA),2)
    CXX = hipcc
    CXXFLAGS = -O3 -std=c++11 -g $(OPENMP) $(CUDACXXFLAGS)
    LD = hipcc -std=c++11
    CC = clang
    CUDACXXFLAGS = -DUSE_CUDA -DAMD_GPU
    CUDALDFLAGS = 
    OPENMP = 
else
    #CMD="g++ -O3 -std=c++11"
    #CMD="g++ -g -O0 -C -std=c++11"
    CXX = icpc
    LD = icpc -std=c++11
    CC = icc
    ifeq ($(MACHINE),amdvega)
        CXX = clang++
        LD = clang++ -std=c++11
        CC = clang
    endif
    ifeq ($(MACHINE),davide)
        CXX = clang++
        LD = clang++ -std=c++11
        CC = clang
    endif
    CXXFLAGS = -O3 -std=c++11 -axCORE-AVX512 -g $(OPENMP) -qopt-report=5 $(CUDACXXFLAGS)
    CUDACXXFLAGS = -DDEBALLOC
    NVCCFLAGS = 
    CUDACXXFLAGS = 
    CUDALDFLAGS =
    CUDAINCLUDE =
endif

ifeq ($(MACHINE),marconi)
    INCLUDEFLAGS = -I/cineca/prod/build/libraries/cgal/4.10/intel--pe-xe-2017--binary/BA_WORK/cgal-4.10/Polyhedron_IO/include/ -I/cineca/prod/opt/libraries/cgal/4.10/intel--pe-xe-2017--binary/include -I/cineca/prod/build/libraries/cgal/4.10/intel--pe-xe-2017--binary/BA_WORK/cgal-4.10/Kernel_23/include/CGAL -I /cineca/prod/build/libraries/cgal/4.10/intel--pe-xe-2017--binary/BA_WORK/cgal-4.10/Polyhedron_IO/include/ -I/cineca/prod/opt/libraries/cgal/4.10/intel--pe-xe-2017--binary/include -I/cineca/prod/opt/libraries/boost/1.61.0/intelmpi--2017--binary/include
    LDFLAGS = -L/cineca/prod/opt/libraries/cgal/4.10/intel--pe-xe-2017--binary/lib64 -lCGAL -L/cineca/prod/opt/libraries/boost/1.61.0/intelmpi--2017--binary/lib -lboost_thread -lboost_system -cxxlib -lboost_chrono -lboost_date_time -lboost_atomic $(OPENMP)
endif
ifeq ($(MACHINE),galileo)
    INCLUDEFLAGS = -I/galileo/home/userinternal/fsalvado/CODE/miniamr/dep/cgal/4.10/cgal-releases-CGAL-4.10/include/CGAL/Polyhedron_IO/include/ -I/galileo/home/userinternal/fsalvado/CODE/miniamr/dep/cgal/4.10/INSTALL_cgal-releases-CGAL-4.10/include/ -I/galileo/home/userinternal/fsalvado/CODE/miniamr/dep/cgal/4.10/cgal-releases-CGAL-4.10/include/CGAL/Kernel_23/include/CGAL -I/cineca/prod/opt/libraries/boost/1.61.0/intelmpi--2017--binary/include $(CUDAINCLUDE)
    LDFLAGS = -L/galileo/home/userinternal/fsalvado/CODE/miniamr/dep/cgal/4.10/INSTALL_cgal-releases-CGAL-4.10/lib64 -lCGAL -L/cineca/prod/opt/libraries/boost/1.61.0/intelmpi--2017--binary/lib -lboost_thread -lboost_system -cxxlib -lboost_chrono -lboost_date_time -lboost_atomic $(OPENMP)
endif
ifeq ($(MACHINE),amdvega)
    INCLUDEFLAGS = -I/home/franz/CODE/miniamr/dep/cgal/cgal-releases-CGAL-4.10.2/Polyhedron_IO/include/ -I/home/franz/CODE/miniamr/dep/cgal/INSTALL_cgal-releases-CGAL-4.10.2/include/ -I./cgal/cgal-releases-CGAL-4.10.2/Kernel_23/include/CGAL 
    LDFLAGS = -L/usr/lib/x86_64-linux-gnu/ -lCGAL -lboost_thread -lboost_system -lboost_chrono -lboost_date_time -lboost_atomic $(OPENMP)
endif
ifeq ($(MACHINE),davide)
    INCLUDEFLAGS = -I/davide/home/userinternal/fsalvado/CODE/miniamr/dep/CGAL-4.10/Polyhedron_IO/include/ -I/davide/home/userinternal/fsalvado/CODE/miniamr/dep/INSTALL_CGAL-4.10/include/ -I./davide/home/userinternal/fsalvado/CODE/miniamr/dep/CGAL-4.10/Kernel_23/include/CGAL -I/davide/home/userinternal/fsalvado/CODE/miniamr/dep/BUILD_boost_1_58_0/include
    LDFLAGS = -L/davide/home/userinternal/fsalvado/CODE/miniamr/dep/INSTALL_CGAL-4.10/lib64 -lCGAL -L/davide/home/userinternal/fsalvado/CODE/miniamr/dep/BUILD_boost_1_58_0/lib -lboost_thread -lboost_system -lboost_chrono -lboost_date_time -lboost_atomic $(OPENMP)
endif

.PHONY: default all clean cleandata

default: $(TARGET)
all: default

SRCS    = $(filter-out $(wildcard *.hipref.cpp), $(wildcard *.cpp))
OBJECTS = $(patsubst %.cpp, %.o, $(SRCS))
ifeq ($(ENABLE_CUDA),1)
    CUDA_OBJECTS = $(patsubst %.cu, %.o, $(wildcard *.cu))
    DEBALLOC = debAlloc.o
else ifeq ($(ENABLE_CUDA),2)
    CUDA_OBJECTS = $(patsubst %.hipref.cpp, %.o, $(wildcard *.hipref.cpp))
    DEBALLOC = 
else
    CUDA_OBJECTS = 
    DEBALLOC = debAlloc.o
endif
HEADERS = $(wildcard *.h)

debAlloc.o: debAlloc.c makefile
	$(CC) -c -DDEBUG_MEMORY debAlloc.c

%.o: %.cpp makefile $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

%.o: %.cu makefile $(HEADERS)
	nvcc $(NVCCFLAGS) $(CUDACXXFLAGS) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

%.o: %.hipref.cpp makefile $(HEADERS)
	hipcc $(NVCCFLAGS) $(CUDACXXFLAGS) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS) $(CUDA_OBJECTS)

$(TARGET): $(OBJECTS) $(CUDA_OBJECTS) $(DEBALLOC)
	@echo "GPU OBJECTS: " $(CUDA_OBJECTS)
	$(LD) $(OBJECTS) $(CUDA_OBJECTS) $(DEBALLOC) $(LDFLAGS) $(INCLUDEFLAGS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)

cleandata:
	-rm -f *.vtr *.vtm output.tar
