ADMB_HOME=/home/mfcl/MFCL/libs/admb/g484-nolp-opt-threaded
PVM_HOME=/home/mfcl/MFCL/libs/pvm3
QD_HOME=/home/mfcl/MFCL/libs/qd-2.3.17
OPENBLAS_HOME=/home/mfcl/MFCL/libs/OpenBLAS-0.2.19

#		*Translator Definitions*
CC= gcc
.SUFFIXES: .cpp .s64
.SUFFIXES: .cpp .s64

TARGET=mfclsdbg64

#     -I/oldhome/dave/src/admb07/admbthreads -I/usr/local/include  \
    
.cpp.s64 : 
	 $(CC) -DNO_MY_DOUBLE_TYPE -fPIC -Wall -fpermissive  -fno-omit-frame-pointer -fPIE  -O0 -femit-class-debug-always  -pthread -o$*.s64   -c -std=gnu++11  \
     -g \
     -DTHREAD_SAFE -Dlinux -DUSE_ADMBTHREAD -DUSE_PTHREADS -D__GNUDOS__ -I. \
     -I$(PVM_HOME)/include -I$(PVM_HOME)/include \
     -I$(QD_HOME)/include \
     -I$(QD_HOME)/include/qd \
     -I$(OPENBLAS_HOME)/include \
     -I$(ADMB_HOME)/include   $<

#		*List Macros*
EXE_dependencies =  \
 indepvars.s64 \
 ils_qr.s64 \
 print_survey.s64 \
 natmort.s64 \
 lapack_inv.s64 \
 recinpop_standard.s64 \
 rationalize_movement.s64 \
 interface.s64 \
 eq2.s64 \
 popes_approx.s64 \
 kludged_equil_surv.s64 \
 nrmissing_effort.s64 \
 setgroup_1.s64 \
 new-len-self-scaling-multiomial.s64 \
 new-wght-self-scaling-multiomial.s64 \
 dfBetai.s64 \
 fvar_a29.s64 \
 ptagfit_ss3.s64 \
 tagfit_ss3.s64 \
 newpredcatch.s64 \
 newway.s64 \
 newmovement2.s64 \
 getcorr1.s64 \
 vsm.s64 \
 newmaturity.s64 \
 convert_matlen_to_matage.s64 \
 sim_realtag_pd.s64 \
 old_new_cross_derivs.s64 \
 oldversion2len.s64 \
 multcoff.s64 \
 makebig2.s64 \
 lognormal_mult_ll_cs.s64 \
 bandedsymstuff.s64 \
 bandedsymlu.s64 \
 get_orth_weights.s64 \
 recinpop_orth.s64 \
 get_orth_poly_info.s64 \
 new_orthogonal_recr.s64 \
 myorth.s64 \
 ssonly.s64 \
 lowtridmatrix.s64 \
 testnewl3.s64 \
 do_all_for_empirical_autocorrelated_bh.s64 \
 sim_tag_pd.s64 \
 sim_xpooltag.s64 \
 get_tag_year.s64 \
 simulation_mode.s64 \
 natural_mortality_spline.s64 \
 minvmult.s64 \
 testsvdinv_noeff.s64 \
 recrpen.s64 \
 df32fun.s64 \
 df22fun.s64 \
 wght_tail_compress.s64 \
 len_tail_compress.s64 \
 learner_code.s64 \
 newfishsel.s64 \
 new_cross_derivs.s64 \
 ssmulttridag.s64 \
 tridiagonal_dmatrix.s64 \
 dmat15.s64 \
 wght_dm_nore.s64 \
 len_dm_nore.s64 \
 version2_len_self_scaling_multinomial_nore.s64 \
 version2_wght_self_scaling_multinomial_nore.s64 \
 version3_len_self_scaling_multinomial_re_multi_rho_multi_var.s64 \
 version3_wght_self_scaling_multinomial_re_multi_rho_multi_var.s64 \
 version2_wght_self_scaling_multinomial_re_multi_rho_multi_var.s64 \
 choleski_check.s64 \
 fvar_m58.s64 \
 lognormal_multinomial4.s64 \
 new_lognormal_multinomial5.s64 \
 phi_newer.s64 \
 minimizer.s64 \
 self_scaling_mult_re.s64 \
 size.s64 \
 newlbsel.s64 \
 new_incident_calc.s64 \
 no_spline.s64 \
 spline.s64 \
 setin.s64 \
 test_squareft.s64 \
 tc_weight_logistic.s64 \
 tc_length_logistic.s64 \
 new_weight_logistic.s64 \
 length_logistic.s64 \
 lwsim.s64  \
 squareft_t.s64  \
 agelength.s64  \
 incident_calc.s64  \
 setcomm3.s64  \
 selbreaks.s64  \
 multspp_tagfit.s64  \
 test_tag_report.s64 \
 newcatcheq.s64 \
 threadtot.s64 \
 xthread_stuff.s64 \
 ythread_stuff.s64 \
 thread_stuff2.s64 \
 thread_stuff3.s64 \
 thread_stuff4.s64 \
 thread_stuff5.s64 \
 thread_stuff6.s64 \
 thread_stuff7.s64 \
 thread_stuff8.s64 \
 threaded_tag3.s64 \
 thread_stuff.s64 \
 fvar_m24.s64 \
 plotstuff.s64 \
 ivec8.s64 \
 gradcforproj.s64 \
 rshort3.s64 \
 newmult.s64 \
 newfmin.s64 \
 newm_io3.s64 \
 fnorm.s64 \
 real_ddnrc3.s64 \
 mfclimplicit.s64 \
 alllengthsel.s64 \
 mfclthread.s64 \
 newrshimp.s64 \
 newrshimp_experiment.s64 \
 fitqimp.s64 \
 nopenalties.s64 \
 goodpen.s64 \
 neworth.s64 \
 testeq1.s64 \
 testeq2.s64 \
 tx.s64 \
 xml.s64 \
 kludge.s64 \
 fishseli.s64 \
 test2.s64 \
 fitcat.s64 \
 nrcatch3x.s64 \
 nrcatch3.s64 \
 nrcatch4.s64 \
 getinp2.s64 \
 avcatfit.s64 \
 avcafit-ms.s64 \
 lesmatrix.s64 \
 fp_rep.s64 \
 yield_bh.s64 \
 yieldbhp.s64 \
 test_msy.s64 \
 testmsy2.s64 \
 flagset.s64 \
 dirmulft.s64 \
 plotN0.s64 \
 onevar.s64 \
 readtag.s64 \
 tagfitsu.s64 \
 tag7.s64 \
 tag8.s64 \
 varcatch.s64 \
 normaliz.s64 \
 normaliz2.s64 \
 setuppvm.s64 \
 equilib.s64 \
 newl6.s64 \
 newl8.s64 \
 xequilib.s64 \
 wfast96.s64 \
 tagtmp.s64 \
 sqrwght.s64 \
 tagfit.s64 \
 ptagfit.s64 \
 grpcatch.s64 \
 biodynam.s64 \
 mlcalc.s64 \
 plot.s64 \
 yukio.s64 \
 optalloc.s64 \
 optfile.s64 \
 kalcat.s64 \
 catbyyr.s64 \
 grptag.s64 \
 pooltag1.s64 \
 xpooltag.s64 \
 diffrout.s64 \
 tag3.s64 \
 www.s64 \
 xwww2.s64 \
 tag4.s64 \
 tag5.s64 \
 rsh3imp.s64 \
 catimp.s64 \
 printfrq.s64 \
 optmatch.s64 \
 lbcutoff.s64 \
 lbselclc.s64 \
 scbound.s64 \
 scsetin1.s64 \
 scsetin.s64 \
 scsetinp.s64 \
 scset1.s64 \
 scsetp.s64 \
 par_ofio.s64 \
 lmul_io4.s64 \
 alb_ctl.s64 \
 nmsort.s64 \
 newcod.s64 \
 derch.s64 \
 readnew.s64 \
 simrep.s64 \
 ests.s64 \
 alldpen2.s64 \
 mfexp.s64 \
 getpath.s64 \
 newgradc.s64 \
 newgradc_noeff.s64 \
 alldevpn.s64 \
 callpen.s64 \
 set2.s64 \
 catlefit.s64 \
 catwt_fit.s64 \
 gradrout.s64 \
 setcomm.s64 \
 setcomm1.s64 \
 setcomm2.s64 \
 setcomm4.s64 \
 setcomm5.s64 \
 setcomm6.s64 \
 setcomm7.s64 \
 exploit.s64 \
 totwtfit.s64 \
 dep_grad.s64 \
 htotcafi.s64 \
 totalfrq.s64 \
 hessrout.s64 \
 robf_fit.s64 \
 rob2_fit.s64 \
 fmeanpen.s64 \
 temppred.s64 \
 veff_dev.s64 \
 rshort1.s64 \
 xrshort3.s64 \
 mwcode.s64 \
 cobbdoug.s64 \
 short2.s64 \
 nfast96.s64 \
 nfast96a.s64 \
 cnfastyr.s64 \
 newl2a.s64 \
 newl2.s64 \
 newl4.s64 \
 newl5.s64 \
 newl9.s64 \
 lmul_io2.s64 \
 lmult.s64 \
 lmult_io.s64 \
 newmaux1.s64 \
 rnaux2.s64 \
 newmaux4.s64 \
 newmaux5.s64 \
 newmau5a.s64 \
 set_ctrl.s64 \
 squareft.s64 \
 vrbioclc.s64 \
 catagfit.s64 \
 crbioclc.s64 \
 envio.s64 \
 admsthread.s64 \
 newgradc_split.s64 \
 mspcode.s64 \
 nnewlan.s64


#		*Explicit Rules*
mfcl: $(EXE_dependencies) 
	$(CC) -fPIE -ggdb3 -pthread  $(EXE_dependencies) -o $(TARGET)  \
  --static -fPIC -std=gnu++11 \
  -L$(QD_HOME)/lib \
  -L$(ADMB_HOME)/lib \
  -L$(PVM_HOME)/lib \
  -L$(OPENBLAS_HOME)/lib \
  -ladt -lads -lstdc++  -ladt -laqds -ladds  -lads -lqd  -lads -lm \
  -ladt -lads -lstdc++  -ladt -laqds -ladds  -lads -lqd  -lads -lm \
  -ladt -lads -lstdc++ \
  -lopenblas  -lgfortran  -lquadmath \
      -ladt \
  -lads -ladds  -lqd  -lads -lm -ladt


#  -L/home/dave/Downloads/pvm3/lib/LINUX64 \

clean:
	rm -f $(EXE_dependencies)

clean-all: clean
	rm -f $(TARGET)
