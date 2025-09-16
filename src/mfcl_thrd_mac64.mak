# Makefile to compile Multifan-CL on Mac
#
#

SHELL :=/bin/bash
ADMB_HOME=../MFCL/lib/admb/g484-nolp-opt-threaded
QD_HOME=../devel/C++/qd-2.3.24
OPENBLAS_HOME:=$(shell brew --prefix openblas)
#FMATH_HOME=$(HOME)/Documents/repository/fmath
GCC_HOME:=$(shell brew --prefix gcc)
CC  :=$(shell ls $(GCC_HOME)/bin |grep ^gcc-[1-9])
CXX :=$(shell ls $(GCC_HOME)/bin |grep ^g++-[1-9])
CPP :=$(shell ls $(GCC_HOME)/bin |grep ^cpp-[1-9]) 
GCCVERSION1:=$(shell $(CC) -dumpversion | cut -f1 -d.)
GCC_MAJOR_VERSION:=$(shell expr $(GCCVERSION1) + 0)

#		*Translator Definitions*
CFLAGS :=  -fpermissive -w
 
CPPFLAGS:= -DTHREAD_SAFE -Dlinux -DUSE_ADMBTHREAD -DUSE_PTHREADS -D__GNUDOS__   -pthread  -DNO_MY_DOUBLE_TYPE -std=gnu++11
ifeq ("$(shell uname)","Darwin")
	CPPFLAGS+= -D_MAC
endif
#LIBS := -lm -lstdc++ -static-libgcc -static-libstdc++ -pthread
LIBS := -lm -lstdc++ -static-libgcc -static-libstdc++ -pthread
ifndef OPT
	OPT=FALSE
endif

ifeq ($(OPT),1)
	OPT=TRUE
endif

ifndef STATIC
	STATIC=TRUE
endif

ifeq ($(STATIC),1)
	STATIC=TRUE
endif

ifeq ($(OPT),TRUE)
	OPTFLAGS :=  -O3   -flto -pipe -fopt-info-vec -fwhole-program
	OPTLIB =	-flto -fwhole-program  -O3
	ifeq ($(STATIC),TRUE)
		OPTFLAGS+=-march=native -mavx2 -mfma
		OPTLIB+=-march=native -mavx2 -mfma
	else
    OPTFLAGS+=-march=sandybridge 
		OPTLIB +=-march=sandybridge 	
	endif
	DEBUGFLAGS =
	CPPFLAGS += -DOPT_LIB
	LIBS += 	-ladt -lado  -laqdo -laddo   
else
	OPTFLAGS =  -O0 -pipe 
	DEBUGFLAGS = -ggdb
	OPTLIB =
	LIBS += 	-ladt -lads -laqds -ladds   
endif

 
#CXXFLAGS = -fvisibility=hidden
CXXFLAGS =
INCLUDES :=	-I/usr/local/include	-I$(QD_HOME)/include 	-I$(QD_HOME)/include/qd 	-I$(OPENBLAS_HOME)/include \
	-I$(ADMB_HOME)/include -I. 
#ifdef FMATH
#	INCLUDES += -I$(FMATH_HOME)
#endif
LDFLAGS :=   -L$(QD_HOME)/lib -L$(ADMB_HOME)/lib -L$(OPENBLAS_HOME)/lib

#### libgomp.a is included YT 2018-07-30
ifeq ($(STATIC),TRUE)
LIBS2 := 	$(OPENBLAS_HOME)/lib/libopenblas.a \
	$(GCC_HOME)/lib/gcc/$(GCC_MAJOR_VERSION)/libgfortran.a \
	$(GCC_HOME)/lib/gcc/$(GCC_MAJOR_VERSION)/libquadmath.a \
	$(GCC_HOME)/lib/gcc/$(GCC_MAJOR_VERSION)/libstdc++.a \
	$(GCC_HOME)/lib/gcc/$(GCC_MAJOR_VERSION)/libgomp.a \
	$(QD_HOME)/lib/libqd.a 
else
	LIBS2 := -lqd   -lopenblas  -lgfortran -lquadmath 
endif
#OPTLIB =	-flto -fwhole-program

ifeq ($(OPT),TRUE)
	O = o64
	ifeq ($(STATIC),TRUE)
		TARGET=mfclo64
	else
		TARGET=mfclo64shlib	
	endif 
else
	O = s64
	ifeq ($(STATIC),TRUE)
		TARGET=mfclsdbg64
	else
		TARGET=mfclsdbgshlib64
	endif
endif
$(warning STATIC=$(STATIC) OPT=$(OPT))
$(warning TARGET=$(TARGET))


.SUFFIXES: .cpp .$(O)

.cpp.$(O) :
	$(CC)  -o$*.$(O)  $(CFLAGS) $(OPTFLAGS) $(DEBUGFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c  $<  \
	
#		*List Macros*
SRC =  \
	indepvars.cpp \
	ils_qr.cpp \
	print_survey.cpp \
	natmort.cpp \
	lapack_inv.cpp \
	recinpop_standard.cpp \
	rationalize_movement.cpp \
	interface.cpp \
	eq2.cpp \
	popes_approx.cpp \
	kludged_equil_surv.cpp \
	nrmissing_effort.cpp \
	setgroup_1.cpp \
	new-len-self-scaling-multiomial.cpp \
	new-wght-self-scaling-multiomial.cpp \
	dfBetai.cpp \
	fvar_a29.cpp \
	ptagfit_ss3.cpp \
	tagfit_ss3.cpp \
	newpredcatch.cpp \
	newway.cpp \
	newmovement2.cpp \
	getcorr1.cpp \
	vsm.cpp \
	newmaturity.cpp \
	convert_matlen_to_matage.cpp \
	sim_realtag_pd.cpp \
	old_new_cross_derivs.cpp \
	oldversion2len.cpp \
	multcoff.cpp \
	makebig2.cpp \
	lognormal_mult_ll_cs.cpp \
	bandedsymstuff.cpp \
	bandedsymlu.cpp \
	get_orth_weights.cpp \
	recinpop_orth.cpp \
	get_orth_poly_info.cpp \
	new_orthogonal_recr.cpp \
	myorth.cpp \
	ssonly.cpp \
	lowtridmatrix.cpp \
	testnewl3.cpp \
	do_all_for_empirical_autocorrelated_bh.cpp \
	sim_tag_pd.cpp \
	sim_xpooltag.cpp \
	get_tag_year.cpp \
	simulation_mode.cpp \
	natural_mortality_spline.cpp \
	minvmult.cpp \
	testsvdinv_noeff.cpp \
	recrpen.cpp \
	df32fun.cpp \
	df22fun.cpp \
	wght_tail_compress.cpp \
	len_tail_compress.cpp \
	learner_code.cpp \
	newfishsel.cpp \
	new_cross_derivs.cpp \
	ssmulttridag.cpp \
	tridiagonal_dmatrix.cpp \
	dmat15.cpp \
	wght_dm_nore.cpp \
	len_dm_nore.cpp \
	version2_len_self_scaling_multinomial_nore.cpp \
	version2_wght_self_scaling_multinomial_nore.cpp \
	version3_len_self_scaling_multinomial_re_multi_rho_multi_var.cpp \
	version3_wght_self_scaling_multinomial_re_multi_rho_multi_var.cpp \
	version2_wght_self_scaling_multinomial_re_multi_rho_multi_var.cpp \
	choleski_check.cpp \
	fvar_m58.cpp \
	lognormal_multinomial4.cpp \
	new_lognormal_multinomial5.cpp \
	phi_newer.cpp \
	minimizer.cpp \
	self_scaling_mult_re.cpp \
	size.cpp \
	newlbsel.cpp \
	new_incident_calc.cpp \
	no_spline.cpp \
	spline.cpp \
	setin.cpp \
	test_squareft.cpp \
	tc_weight_logistic.cpp \
	tc_length_logistic.cpp \
	new_weight_logistic.cpp \
	length_logistic.cpp \
	lwsim.cpp  \
	squareft_t.cpp  \
	agelength.cpp  \
	incident_calc.cpp  \
	setcomm3.cpp  \
	selbreaks.cpp  \
	multspp_tagfit.cpp  \
	test_tag_report.cpp \
	newcatcheq.cpp \
	threadtot.cpp \
	xthread_stuff.cpp \
	ythread_stuff.cpp \
	thread_stuff2.cpp \
	thread_stuff3.cpp \
	thread_stuff4.cpp \
	thread_stuff5.cpp \
	thread_stuff6.cpp \
	thread_stuff7.cpp \
	thread_stuff8.cpp \
	threaded_tag3.cpp \
	thread_stuff.cpp \
	fvar_m24.cpp \
	plotstuff.cpp \
	ivec8.cpp \
	gradcforproj.cpp \
	rshort3.cpp \
	newmult.cpp \
	newfmin.cpp \
	newm_io3.cpp \
	fnorm.cpp \
	real_ddnrc3.cpp \
	mfclimplicit.cpp \
	alllengthsel.cpp \
	mfclthread.cpp \
	newrshimp.cpp \
	newrshimp_experiment.cpp \
	fitqimp.cpp \
	nopenalties.cpp \
	goodpen.cpp \
	neworth.cpp \
	testeq1.cpp \
	testeq2.cpp \
	tx.cpp \
	xml.cpp \
	kludge.cpp \
	fishseli.cpp \
	test2.cpp \
	fitcat.cpp \
	nrcatch3x.cpp \
	nrcatch3.cpp \
	nrcatch4.cpp \
 getinp2.cpp \
 avcatfit.cpp \
 avcafit-ms.cpp \
 lesmatrix.cpp \
 fp_rep.cpp \
 yield_bh.cpp \
 yieldbhp.cpp \
 test_msy.cpp \
 testmsy2.cpp \
 flagset.cpp \
 dirmulft.cpp \
 plotN0.cpp \
 onevar.cpp \
 readtag.cpp \
 tagfitsu.cpp \
 tag7.cpp \
 tag8.cpp \
 varcatch.cpp \
 normaliz.cpp \
 normaliz2.cpp \
 setuppvm.cpp \
 equilib.cpp \
 newl6.cpp \
 newl8.cpp \
 xequilib.cpp \
 wfast96.cpp \
 tagtmp.cpp \
 sqrwght.cpp \
 tagfit.cpp \
 ptagfit.cpp \
 grpcatch.cpp \
 biodynam.cpp \
 mlcalc.cpp \
 plot.cpp \
 yukio.cpp \
 optalloc.cpp \
 optfile.cpp \
 kalcat.cpp \
 catbyyr.cpp \
 grptag.cpp \
 pooltag1.cpp \
 xpooltag.cpp \
 diffrout.cpp \
 tag3.cpp \
 www.cpp \
 xwww2.cpp \
 tag4.cpp \
 tag5.cpp \
 rsh3imp.cpp \
 catimp.cpp \
 printfrq.cpp \
 optmatch.cpp \
 lbcutoff.cpp \
 lbselclc.cpp \
 scbound.cpp \
 scsetin1.cpp \
 scsetin.cpp \
 scsetinp.cpp \
 scset1.cpp \
 scsetp.cpp \
 par_ofio.cpp \
 lmul_io4.cpp \
 alb_ctl.cpp \
 nmsort.cpp \
 newcod.cpp \
 derch.cpp \
 readnew.cpp \
 simrep.cpp \
 ests.cpp \
 alldpen2.cpp \
 mfexp.cpp \
 getpath.cpp \
 newgradc.cpp \
 newgradc_noeff.cpp \
 alldevpn.cpp \
 callpen.cpp \
 set2.cpp \
 catlefit.cpp \
 catwt_fit.cpp \
 gradrout.cpp \
 setcomm.cpp \
 setcomm1.cpp \
 setcomm2.cpp \
 setcomm4.cpp \
 setcomm5.cpp \
	setcomm6.cpp \
	setcomm7.cpp \
 exploit.cpp \
 totwtfit.cpp \
 dep_grad.cpp \
 htotcafi.cpp \
 totalfrq.cpp \
 hessrout.cpp \
 robf_fit.cpp \
 rob2_fit.cpp \
 fmeanpen.cpp \
 temppred.cpp \
 veff_dev.cpp \
 rshort1.cpp \
 xrshort3.cpp \
 mwcode.cpp \
 cobbdoug.cpp \
 short2.cpp \
 nfast96.cpp \
 nfast96a.cpp \
 cnfastyr.cpp \
 newl2a.cpp \
 newl2.cpp \
 newl4.cpp \
 newl5.cpp \
 newl9.cpp \
 lmul_io2.cpp \
 lmult.cpp \
 lmult_io.cpp \
 newmaux1.cpp \
 rnaux2.cpp \
 newmaux4.cpp \
 newmaux5.cpp \
 newmau5a.cpp \
 set_ctrl.cpp \
 squareft.cpp \
 vrbioclc.cpp \
 catagfit.cpp \
 crbioclc.cpp \
 envio.cpp \
 admsthread.cpp \
 newgradc_split.cpp \
 mspcode.cpp \
 nnewlan.cpp

EXE_dependencies =$(SRC:.cpp=.$(O))

dep = $(SRC:.cpp=.d) # one dependency file for each source file

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
%.d: %.cpp
	@$(CPP) $(CFLAGS) $(INCLUDES) $< -MM -MT $(@:.d=.$(O)) >$@

#		*Explicit Rules*
mfcl: $(EXE_dependencies)
	$(CXX)  $(EXE_dependencies) -o $(TARGET)  \
	$(LDFLAGS) $(LIBS) $(LIBS2) $(CXXFLAGS) $(OPTLIB)

-include $(dep) # include all dep files in the makefile

.PHONY: cleandep
cleandep:
	$(RM) $(dep)

.PHONY: clean
clean: cleandep
	$(RM) $(EXE_dependencies) # rm -f $(EXE_dependencies)


.PHONY: clean-all
clean-all: clean 
	$(RM) $(TARGET)   # rm -f $(TARGET)
# DO NOT DELETE
