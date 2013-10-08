all: nvtmain

# Define Fortan compiler
#FC= /opt/intel/composerxe-2011.5.220/bin/ia32/ifort
#FC=/opt/pgi/linux86/10.6/bin/pgfortran  # gfortran
#FC=/usr/local/cuda/bin/nvcc #gfortran
FC=gcc

nvtmain: main.cpp
	$(FC) -o out main.cpp -lm # -deviceemu

#cuda_add.o: cuda_add.cu
#	/usr/local/cuda/bin/nvcc -c cuda_add.cu     #-deviceemu

clean: 
	rm hybsimem main.o
