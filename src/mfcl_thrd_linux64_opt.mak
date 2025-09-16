ADMB_HOME=/home/mfcl/MFCL/libs/admb/g484-nolp-opt-threaded
PVM_HOME=/home/mfcl/MFCL/libs/pvm3
QD_HOME=/home/mfcl/MFCL/libs/qd-2.3.17
OPENBLAS_HOME=/home/mfcl/MFCL/libs/OpenBLAS-0.2.19

#		*Translator Definitions*
CC= gcc
.SUFFIXES: .cpp .o64
.SUFFIXES: .cpp .o64

TARGET=mfclo64

#     -I/oldhome/dave/src/admb07/admbthreads -I/usr/local/include  \
    
.cpp.o64 : 
	 $(CC) --static -w -pthread -o$*.o64  -fpermissive  -c -std=gnu++11  \
     -DOPT_LIB -O2 -DNO_MY_DOUBLE_TYPE \
     -DTHREAD_SAFE -Dlinux -DUSE_ADMBTHREAD -DUSE_PTHREADS -D__GNUDOS__ -I. \
     -I$(PVM_HOME)/include -I$(PVM_HOME)/include \
     -I$(QD_HOME)/include \
     -I$(QD_HOME)/include/qd \
     -I$(OPENBLAS_HOME)/include \
     -I$(ADMB_HOME)/include   $<

 
#		*List Macros*
EXE_dependencies =  \
 indepvars.o64 \
 ils_qr.o64 \
 print_survey.o64 \
 natmort.o64 \
 lapack_inv.o64 \
 recinpop_standard.o64 \
 rationalize_movement.o64 \
 interface.o64 \
 eq2.o64 \
 popes_approx.o64 \
 kludged_equil_surv.o64 \
 nrmissing_effort.o64 \
 setgroup_1.o64 \
 new-len-self-scaling-multiomial.o64 \
 new-wght-self-scaling-multiomial.o64 \
 dfBetai.o64 \
 fvar_a29.o64 \
 ptagfit_ss3.o64 \
 tagfit_ss3.o64 \
 newpredcatch.o64 \
 newway.o64 \
 newmovement2.o64 \
 getcorr1.o64 \
 vsm.o64 \
 newmaturity.o64 \
 convert_matlen_to_matage.o64 \
 sim_realtag_pd.o64 \
 old_new_cross_derivs.o64 \
 oldversion2len.o64 \
 multcoff.o64 \
 makebig2.o64 \
 lognormal_mult_ll_cs.o64 \
 bandedsymstuff.o64 \
 bandedsymlu.o64 \
 get_orth_weights.o64 \
 recinpop_orth.o64 \
 get_orth_poly_info.o64 \
 new_orthogonal_recr.o64 \
 myorth.o64 \
 ssonly.o64 \
 lowtridmatrix.o64 \
 testnewl3.o64 \
 do_all_for_empirical_autocorrelated_bh.o64 \
 sim_tag_pd.o64 \
 sim_xpooltag.o64 \
 get_tag_year.o64 \
 simulation_mode.o64 \
 natural_mortality_spline.o64 \
 minvmult.o64 \
 testsvdinv_noeff.o64 \
 recrpen.o64 \
 df32fun.o64 \
 df22fun.o64 \
 wght_tail_compress.o64 \
 len_tail_compress.o64 \
 learner_code.o64 \
 newfishsel.o64 \
 new_cross_derivs.o64 \
 ssmulttridag.o64 \
 tridiagonal_dmatrix.o64 \
 dmat15.o64 \
 wght_dm_nore.o64 \
 len_dm_nore.o64 \
 version2_len_self_scaling_multinomial_nore.o64 \
 version2_wght_self_scaling_multinomial_nore.o64 \
 version3_len_self_scaling_multinomial_re_multi_rho_multi_var.o64 \
 version3_wght_self_scaling_multinomial_re_multi_rho_multi_var.o64 \
 version2_wght_self_scaling_multinomial_re_multi_rho_multi_var.o64 \
 choleski_check.o64 \
 fvar_m58.o64 \
 lognormal_multinomial4.o64 \
 new_lognormal_multinomial5.o64 \
 phi_newer.o64 \
 minimizer.o64 \
 self_scaling_mult_re.o64 \
 size.o64 \
 newlbsel.o64 \
 new_incident_calc.o64 \
 no_spline.o64 \
 spline.o64 \
 setin.o64 \
 test_squareft.o64 \
 tc_weight_logistic.o64 \
 tc_length_logistic.o64 \
 new_weight_logistic.o64 \
 length_logistic.o64 \
 lwsim.o64  \
 squareft_t.o64  \
 agelength.o64  \
 incident_calc.o64  \
 setcomm3.o64  \
 selbreaks.o64  \
 multspp_tagfit.o64  \
 test_tag_report.o64 \
 newcatcheq.o64 \
 threadtot.o64 \
 xthread_stuff.o64 \
 ythread_stuff.o64 \
 thread_stuff2.o64 \
 thread_stuff3.o64 \
 thread_stuff4.o64 \
 thread_stuff5.o64 \
 thread_stuff6.o64 \
 thread_stuff7.o64 \
 thread_stuff8.o64 \
 threaded_tag3.o64 \
 thread_stuff.o64 \
 fvar_m24.o64 \
 plotstuff.o64 \
 ivec8.o64 \
 gradcforproj.o64 \
 rshort3.o64 \
 newmult.o64 \
 newfmin.o64 \
 newm_io3.o64 \
 fnorm.o64 \
 real_ddnrc3.o64 \
 mfclimplicit.o64 \
 alllengthsel.o64 \
 mfclthread.o64 \
 newrshimp.o64 \
 newrshimp_experiment.o64 \
 fitqimp.o64 \
 nopenalties.o64 \
 goodpen.o64 \
 neworth.o64 \
 testeq1.o64 \
 testeq2.o64 \
 tx.o64 \
 xml.o64 \
 kludge.o64 \
 fishseli.o64 \
 test2.o64 \
 fitcat.o64 \
 nrcatch3x.o64 \
 nrcatch3.o64 \
 nrcatch4.o64 \
 getinp2.o64 \
 avcatfit.o64 \
 avcafit-ms.o64 \
 lesmatrix.o64 \
 fp_rep.o64 \
 yield_bh.o64 \
 yieldbhp.o64 \
 test_msy.o64 \
 testmsy2.o64 \
 flagset.o64 \
 dirmulft.o64 \
 plotN0.o64 \
 onevar.o64 \
 readtag.o64 \
 tagfitsu.o64 \
 tag7.o64 \
 tag8.o64 \
 varcatch.o64 \
 normaliz.o64 \
 normaliz2.o64 \
 setuppvm.o64 \
 equilib.o64 \
 newl6.o64 \
 newl8.o64 \
 newl9.o64 \
 xequilib.o64 \
 wfast96.o64 \
 tagtmp.o64 \
 sqrwght.o64 \
 tagfit.o64 \
 ptagfit.o64 \
 grpcatch.o64 \
 biodynam.o64 \
 mlcalc.o64 \
 plot.o64 \
 yukio.o64 \
 optalloc.o64 \
 optfile.o64 \
 kalcat.o64 \
 catbyyr.o64 \
 grptag.o64 \
 pooltag1.o64 \
 xpooltag.o64 \
 diffrout.o64 \
 tag3.o64 \
 www.o64 \
 xwww2.o64 \
 tag4.o64 \
 tag5.o64 \
 rsh3imp.o64 \
 catimp.o64 \
 printfrq.o64 \
 optmatch.o64 \
 lbcutoff.o64 \
 lbselclc.o64 \
 scbound.o64 \
 scsetin1.o64 \
 scsetin.o64 \
 scsetinp.o64 \
 scset1.o64 \
 scsetp.o64 \
 par_ofio.o64 \
 lmul_io4.o64 \
 alb_ctl.o64 \
 nmsort.o64 \
 newcod.o64 \
 derch.o64 \
 readnew.o64 \
 simrep.o64 \
 ests.o64 \
 alldpen2.o64 \
 mfexp.o64 \
 getpath.o64 \
 newgradc.o64 \
 newgradc_noeff.o64 \
 alldevpn.o64 \
 callpen.o64 \
 set2.o64 \
 catlefit.o64 \
 catwt_fit.o64 \
 gradrout.o64 \
 setcomm.o64 \
 setcomm1.o64 \
 setcomm2.o64 \
 setcomm4.o64 \
 setcomm5.o64 \
 setcomm6.o64 \
 setcomm7.o64 \
 exploit.o64 \
 totwtfit.o64 \
 dep_grad.o64 \
 htotcafi.o64 \
 totalfrq.o64 \
 hessrout.o64 \
 robf_fit.o64 \
 rob2_fit.o64 \
 fmeanpen.o64 \
 temppred.o64 \
 veff_dev.o64 \
 rshort1.o64 \
 xrshort3.o64 \
 mwcode.o64 \
 cobbdoug.o64 \
 short2.o64 \
 nfast96.o64 \
 nfast96a.o64 \
 cnfastyr.o64 \
 newl2a.o64 \
 newl2.o64 \
 newl4.o64 \
 newl5.o64 \
 lmul_io2.o64 \
 lmult.o64 \
 lmult_io.o64 \
 newmaux1.o64 \
 rnaux2.o64 \
 newmaux4.o64 \
 newmaux5.o64 \
 newmau5a.o64 \
 set_ctrl.o64 \
 squareft.o64 \
 vrbioclc.o64 \
 catagfit.o64 \
 crbioclc.o64 \
 envio.o64 \
 admsthread.o64 \
 newgradc_split.o64 \
 mspcode.o64 \
 nnewlan.o64


#		*Explicit Rules*
mfcl: $(EXE_dependencies) 
	$(CC) -pthread  $(EXE_dependencies) -o mfclo64  \
  -static \
  -static -std=gnu++11 \
  -L$(QD_HOME)/lib \
  -L$(ADMB_HOME)/lib \
  -L$(PVM_HOME)/lib \
  -L$(OPENBLAS_HOME)/lib \
  -Xlinker -S \
  -ladt -lado -lstdc++  -ladt -laqdo -laddo  -lado -lqd  -lado -lm \
  -ladt -lado -lstdc++  -ladt -laqdo -laddo  -lado -lqd  -lado -lm \
  -lopenblas  -lgfortran -lquadmath\
  -ladt -lado -lstdc++ \
      -ladt \
  -lado -laddo  -lqd  -lado -lm -ladt


#  -L/home/dave/Downloads/pvm3/lib/LINUX64 \
>

clean:
	rm -f $(EXE_dependencies)

clean-all: clean
	rm -f $(TARGET)
