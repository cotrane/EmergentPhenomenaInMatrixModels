#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputFile(int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, double *lambda,
		char *pathout, int *conf,double *mommult, int *EVstore, int *numrep, int **rep, int *diag, double *prop1, double *prop2, int *compC,
		int *compD, int *compEV, int *compComm, int *compComm2, int *compB)
{
	int i;
	FILE *read;
	char config[1000], *repchar=0;
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../HMC.txt");
	read = fopen(filenameread, "r");
	while(fgets(config, 1000, read)!= NULL)
	{
		//# of loops during simulation
		if(!strcmp(config, "LOOPNUMBER\n"))
		{
			check = fgets(config, 1000, read);
			*lnum = atoi(config);
		}
		//matrix size
		if(!strcmp(config, "SIZE\n"))
		{
			check = fgets(config, 1000, read);
			*size = atoi(config);
		}
		//coupling constant: different def for 3MM and 8MM/15MM!!! see startconf
		if(!strcmp(config, "ALPHATILDE\n"))
		{
			check = fgets(config, 1000, read);
			*atilde = atof(config);
		}
		//stepsize of HMC integration step
		if(!strcmp(config, "EPS\n"))
		{
			check = fgets(config, 1000, read);
			*eps = atof(config);
		}
		//lower bound of acceptance rate of states
		if(!strcmp(config, "LACC\n"))
		{
			check = fgets(config, 1000, read);
			*lacc = atof(config);
		}
		//upper bound of acceptance rate of states
		if(!strcmp(config, "UACC\n"))
		{
			check = fgets(config, 1000, read);
			*uacc = atof(config);
		}
		//loop number where computation of observables starts; supposed loop number where system has thermalized
		if(!strcmp(config, "THERM\n"))
		{
			check = fgets(config, 1000, read);
			*therm = atoi(config);
		}
		//# of integration steps in HMC routine
		if(!strcmp(config, "STEPS\n"))
		{
			check = fgets(config, 1000, read);
			*steps = atoi(config);
		}
		//number of loops between mc steps that are used for computation of observables; if high autocorr useful if >1
		if(!strcmp(config, "STEP\n"))
		{
			check = fgets(config, 1000, read);
			*step = atoi(config);
		}
		//parameter in 8MM/15MM in front of symmetric term
		if(!strcmp(config, "LAMBDA\n"))
		{
			check = fgets(config, 1000, read);
			*lambda = atof(config);
		}
		//path to directory where datafiles are saved
		if(!strcmp(config, "PATHOUT\n"))
		{
			check = fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathout, config);
		}
		//which configuration to start with; in startconf called conf!!!
		if(!strcmp(config, "START\n"))
		{
			check = fgets(config, 1000, read);
			*conf = atoi(config);
		}
		//multiplication factor of entries of gaussian random momenta; increase to increase momenta
		if(!strcmp(config, "MOMMULT\n"))
		{
			check = fgets(config, 1000, read);
			*mommult = atof(config);
		}
		//standard 0; with EVonthefly useless but needs to be zero
		if(!strcmp(config, "EVstore\n"))
		{
			check = fgets(config, 1000, read);
			*EVstore = atoi(config);
		}
		//for 8MM/15MM: if startconfiguration are lie group generators; which size, how many different configurations?
		// #-#-#... # specifies size of one conf in term
		if(!strcmp(config, "REPS\n"))
		{
			check = fgets(config, 1000, read);
			int curr_char = 0;
			while(config[curr_char] != '\n')
			{
				if(config[curr_char] == '-')
					(*numrep)++;
				curr_char++;
			}
			(*numrep)++;
			*rep = (int*) malloc(sizeof(int)*(*numrep));
			memset(*rep, 0, sizeof(int)*(*numrep));
			repchar = strtok(config, "-");
			(*rep)[0] = atoi(repchar);
			for(i=1;i<(*numrep);i++)
			{
				repchar = 0;
				repchar = strtok(NULL, "-");
				(*rep)[i] = atoi(repchar);
			}
		}
		//1 if startconfiguration is SU(2)xU(1); 0 otherwise
		if(!strcmp(config, "DIAG\n"))
		{
			check = fgets(config, 1000, read);
			*diag = atoi(config);
		}
		//?
		if(!strcmp(config, "PROP1\n"))
		{
			check = fgets(config, 1000, read);
			*prop1 = atof(config);
		}
		//?
		if(!strcmp(config, "PROP2\n"))
		{
			check = fgets(config, 1000, read);
			*prop2 = atof(config);
		}
		//include computation of C = gamma_a \otensor X_a ? 0 is no, 1 is yes
		if(!strcmp(config, "compC\n"))
		{
			check = fgets(config, 1000, read);
			*compC = atof(config);
		}
		//include computation of D = gamma_a \otensor [X_a,.] ? 0 is no, 1 is yes
		if(!strcmp(config, "compD\n"))
		{
			check = fgets(config, 1000, read);
			*compD = atof(config);
		}
		//include computation of EV distribution
		if(!strcmp(config, "compEV\n"))
		{
			check = fgets(config, 1000, read);
			*compEV = atof(config);
		}
		//include computation of EV Comm
		if(!strcmp(config, "compComm\n"))
		{
			check = fgets(config, 1000, read);
			*compComm = atof(config);
		}
		//include computation of EV Comm2
		if(!strcmp(config, "compComm2\n"))
		{
			check = fgets(config, 1000, read);
			*compComm2 = atof(config);
		}
		//include computation of EV B
		if(!strcmp(config, "compB\n"))
		{
			check = fgets(config, 1000, read);
			*compB = atof(config);
		}
	}

	fclose(read);
	return 0;
}
