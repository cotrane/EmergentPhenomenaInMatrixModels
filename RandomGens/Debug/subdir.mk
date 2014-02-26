################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../gauss_rand.c \
../mttest.c \
../mtwist.c \
../prng.c \
../randistrs.c \
../rdtest.c 

OBJS += \
./gauss_rand.o \
./mttest.o \
./mtwist.o \
./prng.o \
./randistrs.o \
./rdtest.o 

C_DEPS += \
./gauss_rand.d \
./mttest.d \
./mtwist.d \
./prng.d \
./randistrs.d \
./rdtest.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/home/tkaltenbrunner/CLAPACK-3.2.1/INCLUDE -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


