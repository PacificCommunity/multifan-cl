
TARGET = mfclo64.exe

#Visual Studio home.
MSVC_HOME=C:/Program\ Files\ \(x86\)/Microsoft\ Visual\ Studio/2019/Community/VC/Tools/MSVC/14.29.30133
MSVC_INCLUDE=-I$(MSVC_HOME)/include
MSVC_LIBPATH=/LIBPATH:$(MSVC_HOME)/lib

#Platform SDK home.
MSSDK_HOME=C:/Program\ Files\ \(x86\)/Windows\ Kits/10
MSSDK_INCLUDE=-I$(MSSDK_HOME)/Include
MSSDK_LIBPATH=/LIBPATH:$(MSSDK_HOME)/Lib
MSVCSTUFF1="C:/Program Files (x86)/Windows Kits/10/include/10.0.16299.0/ucrt/"
MSVCSTUFF7="C:/Program Files (x86)/Windows Kits/10/include/10.0.16299.0/um/"
MSVCSTUFF8="C:/Program Files (x86)/Windows Kits/10/include/10.0.16299.0/shared/"

#Pthreads home.
PTHREADS_HOME=../pthreads/pthreads4w-code-v2.10.0-rc/pthreads4w-code-02fecc211d626f28e05ecbb0c10f739bd36d6442
PTHREADS_INCLUDE=-I$(PTHREADS_HOME)
PTHREADS_LIBPATH=/LIBPATH:$(PTHREADS_HOME)

#ADModel home.
ADMODEL_HOME=../my_admb07-apr02-13_dave
ADMODEL_INCLUDE=-I$(ADMODEL_HOME)/include
ADMODEL_LIBPATH=/LIBPATH:$(ADMODEL_HOME)/libs

#OPENBLAS_HOME.
OPENBLAS_HOME=../OpenBLAS/OpenBLAS-v0.2.19-Win64-int32
OPENBLAS_INCLUDE=-I$(OPENBLAS_HOME)/include
OPENBLAS_LIBPATH=/LIBPATH:$(OPENBLAS_HOME)/lib

#QD home.
QD_HOME=../QD/QD-2.3.20/QD-2.3.20
QD_INCLUDE=-I$(QD_HOME)/include
QD_LIBPATH=/LIBPATH:$(QD_HOME)/src/.libs

INCLUDE=$(MSVC_INCLUDE) $(MSSDK_INCLUDE) $(ADMODEL_INCLUDE) $(PTHREADS_INCLUDE) $(QD_INCLUDE) $(OPENBLAS_INCLUDE) -I.

LIBPATH=$(PTHREADS_LIBPATH) $(ADMODEL_LIBPATH) $(QD_LIBPATH) $(OPENBLAS_LIBPATH) $(MSVC_LIBPATH) $(MSSDK_LIBPATH)


LIBS1=$(EXE_dependencies) ado64.lib adto64.lib ddado64.lib qdado64.lib \
 libqd.lib  pthreadVCE2.lib libopenblas.dll.a

#		*Translator Definitions*
#Compiler.
CC=cl
CFLAGS=/Ox /permissive -DOPT_LIB -DNO_MY_DOUBLE_TYPE /MT -DUSE_PTHREADS /c /GF /EHsc -D__MSVC32__=8 -DUSE_LAPLACE -DTHREAD_SAFE -DUSE_ADMBTHREAD -D__WIN32__ /std:c++latest $(INCLUDE)

#Linker.
LD=$(CC)
LDFLAGS=/MACHINE:X64 /nologo  /Fe$(TARGET) $(LIBS1) /link /MANIFEST /MANIFESTFILE:gmult.manifest $(LIBPATH) 

TASM = TASM

TLIB = tlib

TLINK = tlink

.SUFFIXES: .cpp .obm64xo

.cpp.obm64xo : 
	     $(CC) $(CFLAGS) -c -I${MSVCSTUFF1} -I${MSVCSTUFF7} -I${MSVCSTUFF8} /Fo$*.obm64xo $<

#		*List Macros*
EXE_dependencies =  \
 indepvars.obm64xo \
 ils_qr.obm64xo \
 print_survey.obm64xo \
 natmort.obm64xo \
 lapack_inv.obm64xo \
 recinpop_standard.obm64xo \
 rationalize_movement.obm64xo \
 interface.obm64xo \
 eq2.obm64xo \
 popes_approx.obm64xo \
 kludged_equil_surv.obm64xo \
 nrmissing_effort.obm64xo \
 setgroup_1.obm64xo \
 new-len-self-scaling-multiomial.obm64xo \
 new-wght-self-scaling-multiomial.obm64xo \
 dfBetai.obm64xo \
 fvar_a29.obm64xo \
 ptagfit_ss3.obm64xo \
 tagfit_ss3.obm64xo \
 newpredcatch.obm64xo \
 newway.obm64xo \
 newmovement2.obm64xo \
 getcorr1.obm64xo \
 vsm.obm64xo \
 newmaturity.obm64xo \
 convert_matlen_to_matage.obm64xo \
 sim_realtag_pd.obm64xo \
 old_new_cross_derivs.obm64xo \
 oldversion2len.obm64xo \
 multcoff.obm64xo \
 makebig2.obm64xo \
 lognormal_mult_ll_cs.obm64xo \
 bandedsymstuff.obm64xo \
 bandedsymlu.obm64xo \
 get_orth_weights.obm64xo \
 recinpop_orth.obm64xo \
 get_orth_poly_info.obm64xo \
 new_orthogonal_recr.obm64xo \
 myorth.obm64xo \
 ssonly.obm64xo \
 lowtridmatrix.obm64xo \
 testnewl3.obm64xo \
 do_all_for_empirical_autocorrelated_bh.obm64xo \
 sim_tag_pd.obm64xo \
 sim_xpooltag.obm64xo \
 get_tag_year.obm64xo \
 simulation_mode.obm64xo \
 natural_mortality_spline.obm64xo \
 minvmult.obm64xo \
 testsvdinv_noeff.obm64xo \
 recrpen.obm64xo \
 df32fun.obm64xo \
 df22fun.obm64xo \
 wght_tail_compress.obm64xo \
 len_tail_compress.obm64xo \
 learner_code.obm64xo \
 newfishsel.obm64xo \
 new_cross_derivs.obm64xo \
 ssmulttridag.obm64xo \
 tridiagonal_dmatrix.obm64xo \
 dmat15.obm64xo \
 wght_dm_nore.obm64xo \
 len_dm_nore.obm64xo \
 version2_len_self_scaling_multinomial_nore.obm64xo \
 version2_wght_self_scaling_multinomial_nore.obm64xo \
 version3_len_self_scaling_multinomial_re_multi_rho_multi_var.obm64xo \
 version3_wght_self_scaling_multinomial_re_multi_rho_multi_var.obm64xo \
 version2_wght_self_scaling_multinomial_re_multi_rho_multi_var.obm64xo \
 choleski_check.obm64xo \
 fvar_m58.obm64xo \
 lognormal_multinomial4.obm64xo \
 new_lognormal_multinomial5.obm64xo \
 phi_newer.obm64xo \
 minimizer.obm64xo \
 self_scaling_mult_re.obm64xo \
 size.obm64xo \
 newlbsel.obm64xo \
 new_incident_calc.obm64xo \
 no_spline.obm64xo \
 spline.obm64xo \
 setin.obm64xo \
 test_squareft.obm64xo \
 tc_weight_logistic.obm64xo \
 tc_length_logistic.obm64xo \
 new_weight_logistic.obm64xo \
 length_logistic.obm64xo \
 lwsim.obm64xo  \
 squareft_t.obm64xo  \
 agelength.obm64xo  \
 incident_calc.obm64xo  \
 setcomm3.obm64xo  \
 selbreaks.obm64xo  \
 multspp_tagfit.obm64xo  \
 test_tag_report.obm64xo \
 newcatcheq.obm64xo \
 threadtot.obm64xo \
 xthread_stuff.obm64xo \
 ythread_stuff.obm64xo \
 thread_stuff2.obm64xo \
 thread_stuff3.obm64xo \
 thread_stuff4.obm64xo \
 thread_stuff5.obm64xo \
 thread_stuff6.obm64xo \
 thread_stuff7.obm64xo \
 thread_stuff8.obm64xo \
 threaded_tag3.obm64xo \
 thread_stuff.obm64xo \
 fvar_m24.obm64xo \
 plotstuff.obm64xo \
 ivec8.obm64xo \
 gradcforproj.obm64xo \
 rshort3.obm64xo \
 newmult.obm64xo \
 newfmin.obm64xo \
 newm_io3.obm64xo \
 fnorm.obm64xo \
 real_ddnrc3.obm64xo \
 mfclimplicit.obm64xo \
 alllengthsel.obm64xo \
 mfclthread.obm64xo \
 newrshimp.obm64xo \
 newrshimp_experiment.obm64xo \
 fitqimp.obm64xo \
 nopenalties.obm64xo \
 goodpen.obm64xo \
 neworth.obm64xo \
 testeq1.obm64xo \
 testeq2.obm64xo \
 tx.obm64xo \
 xml.obm64xo \
 kludge.obm64xo \
 fishseli.obm64xo \
 test2.obm64xo \
 fitcat.obm64xo \
 nrcatch3x.obm64xo \
 nrcatch3.obm64xo \
 nrcatch4.obm64xo \
 getinp2.obm64xo \
 avcatfit.obm64xo \
 avcafit-ms.obm64xo \
 lesmatrix.obm64xo \
 fp_rep.obm64xo \
 yield_bh.obm64xo \
 yieldbhp.obm64xo \
 test_msy.obm64xo \
 testmsy2.obm64xo \
 flagset.obm64xo \
 dirmulft.obm64xo \
 plotN0.obm64xo \
 onevar.obm64xo \
 readtag.obm64xo \
 tagfitsu.obm64xo \
 tag7.obm64xo \
 tag8.obm64xo \
 varcatch.obm64xo \
 normaliz.obm64xo \
 normaliz2.obm64xo \
 setuppvm.obm64xo \
 equilib.obm64xo \
 newl6.obm64xo \
 newl8.obm64xo \
 xequilib.obm64xo \
 wfast96.obm64xo \
 tagtmp.obm64xo \
 sqrwght.obm64xo \
 tagfit.obm64xo \
 ptagfit.obm64xo \
 grpcatch.obm64xo \
 biodynam.obm64xo \
 mlcalc.obm64xo \
 plot.obm64xo \
 yukio.obm64xo \
 optalloc.obm64xo \
 optfile.obm64xo \
 kalcat.obm64xo \
 catbyyr.obm64xo \
 grptag.obm64xo \
 pooltag1.obm64xo \
 xpooltag.obm64xo \
 diffrout.obm64xo \
 tag3.obm64xo \
 www.obm64xo \
 xwww2.obm64xo \
 tag4.obm64xo \
 tag5.obm64xo \
 rsh3imp.obm64xo \
 catimp.obm64xo \
 printfrq.obm64xo \
 optmatch.obm64xo \
 lbcutoff.obm64xo \
 lbselclc.obm64xo \
 scbound.obm64xo \
 scsetin1.obm64xo \
 scsetin.obm64xo \
 scsetinp.obm64xo \
 scset1.obm64xo \
 scsetp.obm64xo \
 par_ofio.obm64xo \
 lmul_io4.obm64xo \
 alb_ctl.obm64xo \
 nmsort.obm64xo \
 newcod.obm64xo \
 derch.obm64xo \
 readnew.obm64xo \
 simrep.obm64xo \
 ests.obm64xo \
 alldpen2.obm64xo \
 mfexp.obm64xo \
 getpath.obm64xo \
 newgradc.obm64xo \
 newgradc_noeff.obm64xo \
 alldevpn.obm64xo \
 callpen.obm64xo \
 set2.obm64xo \
 catlefit.obm64xo \
 catwt_fit.obm64xo \
 gradrout.obm64xo \
 setcomm.obm64xo \
 setcomm1.obm64xo \
 setcomm2.obm64xo \
 setcomm4.obm64xo \
 setcomm5.obm64xo \
 setcomm6.obm64xo \
 setcomm7.obm64xo \
 exploit.obm64xo \
 totwtfit.obm64xo \
 dep_grad.obm64xo \
 htotcafi.obm64xo \
 totalfrq.obm64xo \
 hessrout.obm64xo \
 robf_fit.obm64xo \
 rob2_fit.obm64xo \
 fmeanpen.obm64xo \
 temppred.obm64xo \
 veff_dev.obm64xo \
 rshort1.obm64xo \
 xrshort3.obm64xo \
 mwcode.obm64xo \
 cobbdoug.obm64xo \
 short2.obm64xo \
 nfast96.obm64xo \
 nfast96a.obm64xo \
 cnfastyr.obm64xo \
 newl2a.obm64xo \
 newl2.obm64xo \
 newl4.obm64xo \
 newl5.obm64xo \
 newl9.obm64xo \
 lmul_io2.obm64xo \
 lmult.obm64xo \
 lmult_io.obm64xo \
 newmaux1.obm64xo \
 rnaux2.obm64xo \
 newmaux4.obm64xo \
 newmaux5.obm64xo \
 newmau5a.obm64xo \
 set_ctrl.obm64xo \
 squareft.obm64xo \
 vrbioclc.obm64xo \
 catagfit.obm64xo \
 crbioclc.obm64xo \
 envio.obm64xo \
 admsthread.obm64xo \
 newgradc_split.obm64xo \
 mspcode.obm64xo \
 nnewlan.obm64xo

#		*Explicit Rules*
all: $(TARGET)

$(TARGET): $(EXE_dependencies) 
	$(LD) $(LDFLAGS)

install: $(TARGET)
	cp $(TARGET) ../bin

clean:
	rm -rf $(EXE_dependencies) $(TARGET)

cleanobj:
	rm -rf $(EXE_dependencies)
