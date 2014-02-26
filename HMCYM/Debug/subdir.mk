################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ActionYM.c \
../Gauss.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../Traces.c \
../computeGdd.c \
../deltaYM.c \
../genD.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c 

OBJS += \
./ActionYM.o \
./Gauss.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./Traces.o \
./computeGdd.o \
./deltaYM.o \
./genD.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o 

C_DEPS += \
./ActionYM.d \
./Gauss.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./Traces.d \
./computeGdd.d \
./deltaYM.d \
./genD.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${CLAPACKdir}/INCLUDE -I${WRKSPACEdir}/RandomGens -I${WRKSPACEdir}/GammaMatrices -I${WRKSPACEdir}/eigenvalue -I${WRKSPACEdir}/MatrixMan -I${WRKSPACEdir}/RepsSUn -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


