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
all: HMCYM_MPI

# Tool invocations
HMCYM_MPI: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	mpicc -L${WRKSPACEdir}/eigenvalue/Release -L${WRKSPACEdir}/RandomGens/Release -L${WRKSPACEdir}/MatrixMan/Release -L${WRKSPACEdir}/RepsSUn/Release -L${WRKSPACEdir}/GammaMatrices/Release -L${CLAPACKdir} -o "HMCYM_MPI" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) HMCYM_MPI
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
