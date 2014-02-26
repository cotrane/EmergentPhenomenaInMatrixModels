################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../GenDiagonal.c \
../GenMultiplets.c \
../GenOffDiagonal.c \
../GenStates.c \
../Main.c \
../Ops.c \
../WeylDimFormula.c 

OBJS += \
./GenDiagonal.o \
./GenMultiplets.o \
./GenOffDiagonal.o \
./GenStates.o \
./Main.o \
./Ops.o \
./WeylDimFormula.o 

C_DEPS += \
./GenDiagonal.d \
./GenMultiplets.d \
./GenOffDiagonal.d \
./GenStates.d \
./Main.d \
./Ops.d \
./WeylDimFormula.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${CLAPACKdir}/INCLUDE -I${WRKSPACEdir}/MatrixMan/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


