################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../3MMwrap.c \
../Action3MM.c \
../Action8MM.c \
../Gauss.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../S1.c \
../action.c \
../deltaS1.c \
../deltaX3MM.c \
../deltaX8MM.c \
../findStructConst.c \
../genD.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c \
../startconf.c \
../symmstructconst.c \
../thecode.c 

OBJS += \
./3MMwrap.o \
./Action3MM.o \
./Action8MM.o \
./Gauss.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./S1.o \
./action.o \
./deltaS1.o \
./deltaX3MM.o \
./deltaX8MM.o \
./findStructConst.o \
./genD.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o \
./startconf.o \
./symmstructconst.o \
./thecode.o 

C_DEPS += \
./3MMwrap.d \
./Action3MM.d \
./Action8MM.d \
./Gauss.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./S1.d \
./action.d \
./deltaS1.d \
./deltaX3MM.d \
./deltaX8MM.d \
./findStructConst.d \
./genD.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d \
./startconf.d \
./symmstructconst.d \
./thecode.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${WRKSPACEdir}/MatrixMan -I${WRKSPACEdir}/RandomGens/ -I${WRKSPACEdir}/GammaMatrices -I${WRKSPACEdir}/eigenvalue -I${WRKSPACEdir}/RepsSUn -I${CLAPACKdir}/INCLUDE -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


