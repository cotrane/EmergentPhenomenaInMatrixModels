################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../AddMat.c \
../AntiComm.c \
../Comm.c \
../Commend.c \
../Inverse.c \
../Multi.c \
../Multilambda.c \
../diagMulti.c \
../printmat.c 

OBJS += \
./AddMat.o \
./AntiComm.o \
./Comm.o \
./Commend.o \
./Inverse.o \
./Multi.o \
./Multilambda.o \
./diagMulti.o \
./printmat.o 

C_DEPS += \
./AddMat.d \
./AntiComm.d \
./Comm.d \
./Commend.d \
./Inverse.d \
./Multi.d \
./Multilambda.d \
./diagMulti.d \
./printmat.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"${CLAPACKdir}/INCLUDE" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


