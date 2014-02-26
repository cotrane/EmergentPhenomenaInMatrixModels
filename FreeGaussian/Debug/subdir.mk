################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Gauss.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../genGauss_NonHermitianMatrix.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c 

OBJS += \
./Gauss.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./genGauss_NonHermitianMatrix.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o 

C_DEPS += \
./Gauss.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./genGauss_NonHermitianMatrix.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${WRKSPACEdir}/RandomGens/ -I${CLAPACKdir}/INCLUDE -I${WRKSPACEdir}/RepsSUn/ -I${WRKSPACEdir}/GammaMatrices/ -I${WRKSPACEdir}/MatrixMan/ -I${WRKSPACEdir}/eigenvalue/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


