################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../1MMwrap.c \
../3MMwrap.c \
../Action1MM.c \
../Action3MM.c \
../Action8MM.c \
../Gauss.c \
../Mass.c \
../ReadInputCommandLine.c \
../ReadInputFile.c \
../S1.c \
../computesc.c \
../deltaS1.c \
../deltaX3MM.c \
../deltaX8MM.c \
../dtensor.c \
../generate_gaussp.c \
../generate_randphi.c \
../main.c \
../startconf.c \
../structconst.c \
../symmstructconst.c \
../thecode.c \
../timerev.c 

OBJS += \
./1MMwrap.o \
./3MMwrap.o \
./Action1MM.o \
./Action3MM.o \
./Action8MM.o \
./Gauss.o \
./Mass.o \
./ReadInputCommandLine.o \
./ReadInputFile.o \
./S1.o \
./computesc.o \
./deltaS1.o \
./deltaX3MM.o \
./deltaX8MM.o \
./dtensor.o \
./generate_gaussp.o \
./generate_randphi.o \
./main.o \
./startconf.o \
./structconst.o \
./symmstructconst.o \
./thecode.o \
./timerev.o 

C_DEPS += \
./1MMwrap.d \
./3MMwrap.d \
./Action1MM.d \
./Action3MM.d \
./Action8MM.d \
./Gauss.d \
./Mass.d \
./ReadInputCommandLine.d \
./ReadInputFile.d \
./S1.d \
./computesc.d \
./deltaS1.d \
./deltaX3MM.d \
./deltaX8MM.d \
./dtensor.d \
./generate_gaussp.d \
./generate_randphi.d \
./main.d \
./startconf.d \
./structconst.d \
./symmstructconst.d \
./thecode.d \
./timerev.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I${WRKSPACEdir}/eigenvalue -I${WRKSPACEdir}/RandomGens -I${WRKSPACEdir}/MatrixMan -I${WRKSPACEdir}/RepsSUn -I${CLAPACKdir}/INCLUDE -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


