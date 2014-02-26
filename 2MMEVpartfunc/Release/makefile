################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: 2MMEVpartfunc

# Tool invocations
2MMEVpartfunc: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -L${CLAPACKdir}/ -L"${WRKSPACEdir}/RandomGens/Release" -L"${WRKSPACEdir}/RepsSUn/Release" -L"${WRKSPACEdir}/GammaMatrices/Release" -L"${WRKSPACEdir}/eigenvalue/Release" -L${WRKSPACEdir}/MatrixMan/Release -o"2MMEVpartfunc" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) 2MMEVpartfunc
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
