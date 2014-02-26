import os
import subprocess
import math

g=0
for i in range(0,1):
	loopnr = 64000
	matsize = 10
	nummat = 2
	g = 0
	mass = 0.010 + 0.001*i
	eps = 0.03
	lacc = 70
	uacc = 90
	therm = 32000
	intsteps = 2
	step = 2
	pathout = "./"
	mult = 1.
	EVstore = 0
	compC = 0
	compD = 1
	compEV = 0
	compComm = 0
	compTraces = 0
	compPhi = 0
	compOffDiag = 0
	compEV = 0

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
	arg22 = "MULT"
	arg24 = "EVstore"
	arg26 = "compC"
	arg28 = "compD"
	arg30 = "compEV"
	arg32 = "compComm"
	arg34 = "compTraces"
	arg36 = "compPhi"
	arg38 = "compAComm"
	arg40 = "NUMMAT"
	arg42 = "compOffDiag"
	arg44 = "compEV"
	if g!=0:
		arg20 = "ParG"
	if mass!=0:
		arg20 = "MASS"

	arg3 = "%d" % loopnr
	arg5 = "%d" % matsize
	arg7 = "%d" % therm
	arg9 = "%d" % step
	arg11 = "%d" % intsteps
	arg13 = "%f" % lacc
	arg15 = "%f" % uacc
	arg17 = "%f" % eps
	arg19 = "%s" % pathout
	arg25 = "%f" % EVstore
	arg27 = "%d" % compC
	arg29 = "%d" % compD
	arg31 = "%d" % compEV
	arg33 = "%d" % compComm
	arg35 = "%d" % compTraces
	arg37 = "%d" % compPhi
	arg39 = "%d" % compAComm
	arg41 = "%d" % nummat
	arg43 = "%d" % compOffDiag
	arg45 = "%d" % compEV
	if g!=0:
		arg21 = "%f" % g
	if mass!=0:
		arg21 = "%f" % mass
	arg23 = "%f" % mult

	print arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, \
		arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, \
		arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg39, \
		arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45

	retcode = subprocess.call(["/home/tkaltenbrunner/workspace/2MMeff", arg1, arg2, arg3, arg1, arg4, arg5, arg1, arg6, arg7, arg1, arg8, arg9, arg1, arg10, arg11, arg1, arg12, arg13, arg1, arg14, arg15, arg1, arg16, arg17, arg1, arg18, arg19, arg1, arg20, arg21, arg1, arg22, arg23, arg1, arg24, arg25, arg1, arg26, arg27, arg1, arg28, arg29, arg1, arg30, arg31, arg1, arg32, arg33, arg1, arg34, arg35, arg1, arg36, arg37, arg1, arg38, arg39, arg1, arg40, arg41, arg1, arg42, arg43, arg1, arg44, arg45])

