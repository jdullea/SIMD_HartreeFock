################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BasisFunction.cpp \
../src/VectorizedHF.cpp 

OBJS += \
./src/BasisFunction.o \
./src/VectorizedHF.o 

CPP_DEPS += \
./src/BasisFunction.d \
./src/VectorizedHF.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/Users/jon/software/install/libint2/include/libint2 -I/Users/jon/software/install/libint2/include -I/opt/local/include/eigen3 -O3 -Wall -c -fmessage-length=0 -std=c++11 -O3 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


