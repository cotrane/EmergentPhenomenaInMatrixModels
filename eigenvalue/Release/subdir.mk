################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../DataDistribution.c \
../EVonthefly.c \
../EVtoDist_Modulus_nH.c \
../EVtoDist_nonHermitian.c \
../EigenV.c \
../FileClose.c \
../FileOpen.c \
../SaveEV.c \
../SaveEV_1.c \
../SaveEVendofflight.c \
../StoreData.c \
../StoreEV.c 

OBJS += \
./DataDistribution.o \
./EVonthefly.o \
./EVtoDist_Modulus_nH.o \
./EVtoDist_nonHermitian.o \
./EigenV.o \
./FileClose.o \
./FileOpen.o \
./SaveEV.o \
./SaveEV_1.o \
./SaveEVendofflight.o \
./StoreData.o \
./StoreEV.o 

C_DEPS += \
./DataDistribution.d \
./EVonthefly.d \
./EVtoDist_Modulus_nH.d \
./EVtoDist_nonHermitian.d \
./EigenV.d \
./FileClose.d \
./FileOpen.d \
./SaveEV.d \
./SaveEV_1.d \
./SaveEVendofflight.d \
./StoreData.d \
./StoreEV.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${CLAPACKdir}/INCLUDE -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


