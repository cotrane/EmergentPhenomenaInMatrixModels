import os
import subprocess
import math

#f = open('/home/tkaltenbrunner/workspace/pars15a3.3.txt', 'r')
for i in range(0,1):
	loopnr = 128000
	matsize = 6
	alpha = 2.10019 + 0.2*i
	eps = 0.03
	lacc = 70
	uacc = 90
	therm = 32000
	intsteps = 2
	step = 1
	lammbda = 0.0 #+ 0.01*i
	start = 0
	evstore = 0
	mommult=1.0
	diag = 0
	reps = 0
	compC = 0
	compD = 0
	compEV = 1
	compComm = 1
	compB = 0
	compEV = 0
	prop1 = 1.0
	prop2 = 1.0
	pathout = "./" % (matsize,matsize)

	arg1 = "-par"
	arg2 = "ALPHATILDE"
	arg4 = "LOOPNUMBER"
	arg6 = "SIZE"
	arg8 = "THERM"
	arg10 = "STEP"
	arg12 = "STEPS"
	arg14 = "LACC"
	arg16 =  "UACC"
	arg18 = "LAMBDA"
	arg20 = "EPS"
	arg22 = "PATHOUT"
	arg24 = "START"
	arg26 = "EVstore"
	arg28 = "MOMMULT"
	arg30 = "REPS"
	arg32 = "DIAG"
	arg34 = "PROP1"
	arg36 = "PROP2"
	arg38 = "compC"
	arg40 = "compD"
	arg42 = "compEV"
	arg44 = "compComm"
	arg46 = "compB"
	arg48 = "compEV"
	arg50 = "prop1"
	arg52 = "prop2"

	arg3 = "%f" % alpha
	arg5 = "%d" % loopnr
	arg7 = "%d" % matsize
	arg9 = "%d" % therm
	arg11 = "%d" % step
	arg13 = "%d" % intsteps
	arg15 = "%f" % lacc
	arg17 = "%f" % uacc
	arg19 = "%f" % lammbda
	arg21 = "%f" % eps
	arg23 = "%s" % pathout
	arg25 = "%d" % start
	arg27 = "%d" % evstore
	arg29 = "%f" % mommult
	arg31 = "%d" % reps
	arg39 = "%d" % compC
	arg41 = "%d" % compD
	arg43 = "%d" % compEV
	arg45 = "%d" % compComm
	arg47 = "%d" % compB
	arg49 = "%d" % compEV
	arg51 = "%d" % prop1
	arg53 = "%d" % prop2

	if(start==2 and diag==1):
		reps = int(f.readline())
		prop1 = float(f.readline())
		prop2 = float(f.readline())
		reps = reps-1
		arg31 = "%s" % reps
		arg33 = "%d" % diag
		arg35 = "%f" % prop1
		arg37 = "%f" % prop2

	if(start==2 and diag==0 and lammbda!=0):
		prop1 = 1#float(f.readline())
		reps = matsize-1
		arg31 = "%s" % reps
		arg33 = "%d" % diag
		arg35 = "%f" % prop1

	if(start==2 and diag==0 and lammbda==0):
		reps = matsize-1
		arg31 = "%s" % reps

	if(start==2 and diag==1):
		print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
			arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27,\
	 		arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg39, arg1, arg40, arg41, \
			arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53

		retcode = subprocess.call(["./HMC", arg1, arg2, arg3, arg1, arg4, arg5,\
	 		arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19,\
	  		arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, \
	  		arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, \
			arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53])

	elif(start==2 and diag==0 and lammbda!=0):
		print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
			arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27,\
	 		arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, \
			arg1, arg44, arg45, arg1, arg46, arg47, arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53

		retcode = subprocess.call(["./HMC", arg1, arg2, arg3, arg1, arg4, arg5,\
	  		arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19,\
	  		arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, \
	  		arg1, arg34, arg35, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, arg1, arg48, arg49, \
			arg1, arg50, arg51, arg1, arg52, arg53])

	elif(start==2 and diag==0 and lammbda==0):
		print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
			arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27,\
	 		arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, \
			arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53

		retcode = subprocess.call(["./HMC", arg1, arg2, arg3, arg1, arg4, arg5,\
	  		arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19,\
	  		arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg38, arg39, \
			arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53])

	else:
		print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
			arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27,\
	 		arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, \
			arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53

		retcode = subprocess.call(["./HMC", arg1, arg2, arg3, arg1, arg4, arg5,\
	  		arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19,\
	  		arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg38, arg39, \
			arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45, arg1, arg46, arg47, arg1, arg48, arg49, arg1, arg50, arg51, arg1, arg52, arg53])

#f.close()


#"mpirun", "-n", "8", "./HMC_MPI"

