################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ScanMatrix.c \
../gammasize.c \
../main.c 

OBJS += \
./ScanMatrix.o \
./gammasize.o \
./main.o 

C_DEPS += \
./ScanMatrix.d \
./gammasize.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${CLAPACKdir}/INCLUDE -I${WRKSPACEdir}/MatrixMan/ -I${WRKSPACEdir}/RepsSUn/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


