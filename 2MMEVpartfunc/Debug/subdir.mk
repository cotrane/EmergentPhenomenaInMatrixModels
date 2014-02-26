################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ActionYM.c \
../ConstructAdjointAction.c \
../Gauss.c \
../GenGaussMatrix.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../Traces.c \
../deltaYM.c \
../findStructConst.c \
../genD.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c 

OBJS += \
./ActionYM.o \
./ConstructAdjointAction.o \
./Gauss.o \
./GenGaussMatrix.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./Traces.o \
./deltaYM.o \
./findStructConst.o \
./genD.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o 

C_DEPS += \
./ActionYM.d \
./ConstructAdjointAction.d \
./Gauss.d \
./GenGaussMatrix.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./Traces.d \
./deltaYM.d \
./findStructConst.d \
./genD.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"/home/tkaltenbrunner/CodeForBadis/MatrixMan" -I"/home/tkaltenbrunner/CodeForBadis/RandomGens" -I"/home/tkaltenbrunner/CodeForBadis/GammaMatrices" -I"/home/tkaltenbrunner/CodeForBadis/eigenvalue" -I"/home/tkaltenbrunner/CodeForBadis/RepsSUn" -I/home/tkaltenbrunner/CLAPACK-3.2.1/INCLUDE -O0 -g3 -Wall -c -fmessage-length=0 -Wno-unused-variable -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


