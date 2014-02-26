import os
import subprocess
import math

for i in range(0,1):
	nummat = 3 + 1*i;
	loopnr =64000
	matsize = 6# + i
	mass = 0.0 #+ 0.01*i
	eps = 0.005
	lacc = 60
	uacc = 80
	therm = 32000
	intsteps = 8
	step=1
	evstore = 0
	compC = 0
	compD = 0
	compEV = 1
	compComm = 0
	compTraces = 0
	compGdd = 0
	pathout = "./"
	liegroup = 2
	start = 1

	arg1 = "-par"
	arg2 = "LOOPNUMBER"
	arg4 = "SIZE"
	arg6 = "THERM"
	arg8 = "STEP"
	arg10 = "STEPS"
	arg12 = "LACC"
	arg14 =  "UACC"
	arg16 = "EPS"
	arg18 = "PATHOUT"
	arg20 = "EVstore"
	arg22 = "NUMMAT"
	arg24 = "LIEGROUP"
	arg26 = "MASS"
	arg28 = "compC"
	arg30 = "compD"
	arg32 = "compEV"
	arg34 = "compComm"
	arg36 = "compTraces"
	arg38 = "START"
	arg40 = "compGdd"

	arg3 = "%d" % loopnr
	arg5 = "%d" % matsize
	arg7 = "%d" % therm
	arg9 = "%d" % step
	arg11 = "%d" % intsteps
	arg13 = "%f" % lacc
	arg15 = "%f" % uacc
	arg17 = "%f" % eps
	arg19 = "%s" % pathout
	arg21 = "%d" % evstore
	arg23 = "%d" % nummat
	arg25 = "%d" % liegroup
	arg27 = "%f" % mass
	arg29 = "%d" % compC
	arg31 = "%d" % compD
	arg33 = "%d" % compEV
	arg35 = "%d" % compComm
	arg37 = "%d" % compTraces
	arg39 = "%d" % start
	arg41 = "%d" % compGdd

	print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
		arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26,	arg27, \
		arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg40, arg1, arg39, arg41

	retcode = subprocess.call(["./HMCYM", arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg40, arg1, arg39, arg41])

