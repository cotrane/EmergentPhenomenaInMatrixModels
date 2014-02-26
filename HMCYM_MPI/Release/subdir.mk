################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ActionYM.c \
../Gauss.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../computeGdd.c \
../deltaYM.c \
../genD.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c \
../printmat.c 

OBJS += \
./ActionYM.o \
./Gauss.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./computeGdd.o \
./deltaYM.o \
./genD.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o \
./printmat.o 

C_DEPS += \
./ActionYM.d \
./Gauss.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./computeGdd.d \
./deltaYM.d \
./genD.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d \
./printmat.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	mpicc -I${CLAPACKdir}/INCLUDE -I${WRKSPACEdir}/RandomGens -I${WRKSPACEdir}/eigenvalue -I${WRKSPACEdir}/RepsSUn -I${WRKSPACEdir}/MatrixMan -I${WRKSPACEdir}/GammaMatrices -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


