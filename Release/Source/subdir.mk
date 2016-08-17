################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Source/HamiltonianBlockSolver.cpp \
../Source/main.cpp \
../Source/szBasis.cpp \
../Source/szHamiltonian.cpp \
../Source/szTransBasis.cpp \
../Source/szTransHamiltonian.cpp \
../Source/utils.cpp 

OBJS += \
./Source/HamiltonianBlockSolver.o \
./Source/main.o \
./Source/szBasis.o \
./Source/szHamiltonian.o \
./Source/szTransBasis.o \
./Source/szTransHamiltonian.o \
./Source/utils.o 

CPP_DEPS += \
./Source/HamiltonianBlockSolver.d \
./Source/main.d \
./Source/szBasis.d \
./Source/szHamiltonian.d \
./Source/szTransBasis.d \
./Source/szTransHamiltonian.d \
./Source/utils.d 


# Each subdirectory must supply rules for building sources it contributes
Source/%.o: ../Source/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


