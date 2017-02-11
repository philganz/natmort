//===================================================================================================
DATA_SECTION
//===================================================================================================
  !!CLASS ofstream mcmc_results("mcmc_results.csv",ios::trunc);

// Read non-simulated data

  // Number of years
  init_int      nyrs
  // Number of ages
  init_int      nages
  // Spawn month
  init_number spawn_fract
  !! spawn_fract = (spawn_fract-1)/12;
  // Recruitment age
  init_int recage
  // Weight-at-age
  init_vector   wt_F(1,nages);
  init_vector   wt_M(1,nages);
  vector        wt(1,nages);
  !! wt = (wt_F+wt_M)/2;  
  // Maturity
  init_vector   mat(1,nages);
  vector        wt_mat(1,nages);
  !! wt_mat = elem_prod(wt_F,mat)/2;
  // Natural mortality covariate
  init_vector obs_m_cov(1,nyrs);

// Read simulated data

  // Longline survey biomass CV
  init_number   obs_srv_biom_CV
  // Longline survey biomass
  init_vector   obs_srv_biom(1,nyrs);
  // Fishery biomass CV
  init_number   obs_fish_biom_CV
  // Fishery biomass
  init_vector   obs_fish_biom(1,nyrs);
  // Longline survey age comps sample size  
  init_number   nsamples_srv_age
  // Longline survey age comps
  init_matrix   obs_ac_srv(1,nyrs,1,nages);
  // Fishery age comps sample size  
  init_number   nsamples_fish_age
  // Fishery age comps
  init_matrix   obs_ac_fish(1,nyrs,1,nages);

  //EOF Marker
  init_int eof;


  !! ad_comm::change_datafile_name("tem.ctl");
// Read parameter phases
  init_int    ph_logR
  init_int    ph_Rdevs
  init_int    ph_Idevs
  init_int    ph_M_0
  init_int    ph_M_1
  init_int    ph_a  
  init_int    ph_B
  init_int    ph_avgF
  init_int    ph_Fdevs
  init_int    ph_q
  init_int    ph_Ssel
  init_int    ph_Fsel
  init_int    ph_phi
  init_int    ph_Mdevs
  init_int    ph_sig
// Read in natural mortality estimation case
  init_int    M_case

// Initialize some counting variables
  int i
  int j
  int header
  !!header = 1;

 LOCAL_CALCS

  //End of file marker
   if(eof==42) cout<<"The data has been read correctly!";
   else { cout <<"What we have here is failure to communicate"<<endl;exit(1);  }

 END_CALCS

//===================================================================================================
PARAMETER_SECTION
//===================================================================================================

//==== Define parameters

// Recruitment/initial abundance parameters
  init_bounded_number       logR(-1,3,ph_logR);
  init_bounded_dev_vector   rec_devs(1,nyrs,-15,15,ph_Rdevs);
  init_bounded_vector       init_devs(2,nages,-5,5,ph_Idevs);

// Fishing mortality
  init_number               log_avg_F(ph_avgF);
  init_bounded_dev_vector   F_devs(1,nyrs,-15,15,ph_Fdevs);

// Natural mortality
  init_bounded_number             log_M_0(-5,0,ph_M_0);
  number                          M_0;
  init_number             log_M_1(ph_M_1);        // For random walk scenarios (log_M_1 = log(M(1))) 
  init_bounded_number     log_phi(-10000,0,ph_phi);
  number                  phi;
  init_number             alpha(ph_a);
  init_number             Beta(ph_B);
  
// Survey catchability
  init_number         log_q_srv(ph_q);

// Survey selectivity
  init_number         a50_srv(ph_Ssel);
  init_number         delta_srv(ph_Ssel);

// Fishery selectivity
  init_number         a50_fish(ph_Fsel);
  init_number         delta_fish(ph_Fsel);

//==== Define dependent numbers/vectors/matrices

// Survey selectivity
  vector  srv_sel(1,nages);

// Fishery selectivity
  vector  fish_sel(1,nages);

// Numbers at age
  matrix        natage(1,nyrs,1,nages);

// Mortality and survival
  vector        Fmort(1,nyrs);
  matrix        F(1,nyrs,1,nages);
  vector        M(1,nyrs);
  matrix        Z(1,nyrs,1,nages);
  matrix        S(1,nyrs,1,nages);
  vector        pred_rec(1,nyrs);

// Predicted values
  matrix        n_srv(1,nyrs,1,nages);
  sdreport_vector        pred_srv(1,nyrs);
  matrix        eac_srv(1,nyrs,1,nages);

  matrix        catage(1,nyrs,1,nages);
  sdreport_vector        pred_catch(1,nyrs);
  matrix        eac_fish(1,nyrs,1,nages);

// Standard deviation estimates
  sdreport_vector       tot_biom(1,nyrs);
  sdreport_vector       spawn_biom(1,nyrs);

// Random effects (and associated sigma)
  init_bounded_number           log_sigma_M(-20,20,ph_sig);
  number                sigma_M;

// Natural mortality as fixed effects vector
//  init_vector   M_devs(1,nyrs,ph_Mdevs);

// Natural mortality deviations as random effects
  random_effects_vector M_devs(1,nyrs,ph_Mdevs);

// Change vector length for walks and correlated scenarios
//  init_vector  M_devs(2,nyrs,ph_Mdevs);
//  random_effects_vector  M_devs(2,nyrs,ph_Mdevs);


// Likelihoods and penalty functions
  number         srv_like;
  number         catch_like;
  number         srv_age_like;
  number         fish_age_like;
  objective_function_value obj_fun;

//===================================================================================================
PROCEDURE_SECTION
//===================================================================================================

  Get_Selectivity();
  Get_Mortality_Rates();
  Get_Numbers_At_Age();
  Get_Predicted_Values();
  Get_Age_Comp();
  Evaluate_Objective_Function();
  if(last_phase()){
  write_base_results();}
  if(mceval_phase()){
  write_mcmc_results();}

//===================================================================================================
FUNCTION Get_Selectivity
//===================================================================================================

//  Survey selectivity
  for (j=1;j<=nages;j++) {
    srv_sel(j)=1/(1+mfexp(-delta_srv*((j)+1-a50_srv)));}

//  Fishery selectivity
  for (j=1;j<=nages;j++) {
    fish_sel(j)=1/(1+mfexp(-delta_fish*((j)+1-a50_fish)));}

//===================================================================================================
FUNCTION Get_Mortality_Rates
//===================================================================================================

// Transformations
  M_0 = mfexp(log_M_0);
  sigma_M = mfexp(log_sigma_M);
  phi = mfexp(log_phi);
 
// Natural mortality  
  // Covariate case
  // Turn beta phase off for estimation without covariate
  // Turn beta and log_M_0 phase off to fix natural mortality
  if(M_case==1){
  for (i=1;i<=nyrs;i++){
  M(i) = M_0+Beta*obs_m_cov(i);}}

  // Uncorrelated deviations
  if(M_case==2){
  for (i=1;i<=nyrs;i++){
  M(i) = M_0+sigma_M*M_devs(i);}}

  // Random walk 
  if(M_case==3){
  M(1)=mfexp(log_M_1);
  for (i=2;i<=nyrs;i++){
  M(i) = M(i-1)+alpha+sigma_M*M_devs(i);}
  M_0 =mean(M);}
 
  // Correlated walk
  if(M_case==4){
  M(1)=mfexp(log_M_1);
  for (i=2;i<=nyrs;i++){
  M(i) = phi*M(i-1)+alpha+sigma_M*M_devs(i);}
  M_0 =mean(M);}

// Fishing mortality
  Fmort = mfexp(log_avg_F + F_devs);
  for (i=1;i<=nyrs;i++){
   F(i) = Fmort(i)*fish_sel;}  

// Total mortality  
  for (i=1;i<=nyrs;i++){
  for (j=1;j<=nages;j++){
   Z(i,j) = F(i,j) + M(i);}}

// Survival  
  for (i=1;i<=nyrs;i++){
  for (j=1;j<=nages;j++){
  S(i,j) = mfexp(-1.0*Z(i,j));}}
  
//===================================================================================================
FUNCTION Get_Numbers_At_Age  
//===================================================================================================

//==== Calculate Numbers at age

// Initial abundance
    for (j=2;j<=nages;j++){
      natage(1,j) = mfexp(init_devs(j));}

// Recruitment
  for (i=1;i<=nyrs;i++){
    natage(i,1) = mfexp(logR+rec_devs(i));}
  
// Abundance after first age and year
  for (i=2;i<=nyrs;i++){
   for (j=2;j<=nages-1;j++){
    natage(i,j) = natage(i-1,j-1)*S(i-1,j-1);}
   natage(i,nages) = (natage(i-1,nages-1)*S(i-1,nages-1))+(natage(i-1,nages)*S(i-1,nages));
  }

//==== Calculate recruitment and total/spawning biomass

  for (i=1;i<=nyrs;i++){
    tot_biom(i) = wt * natage(i);
    spawn_biom(i) = wt_mat * natage(i);
    pred_rec(i) = natage(i,1);}

//===================================================================================================
FUNCTION Get_Predicted_Values  
//===================================================================================================
  
  // Survey 
  // Numbers-at-age
  for (i=1;i<=nyrs;i++){
   for (j=1;j<=nages;j++){
    n_srv(i,j) = mfexp(log_q_srv)*srv_sel(j)*natage(i,j);}}
  
  // Biomass
  for (i=1;i<=nyrs;i++){
    pred_srv(i) = n_srv(i)*wt;}

  // Fishery 
  // Numbers-at-age
  for (i=1;i<=nyrs;i++){
   catage(i) = elem_div(elem_prod(elem_prod(natage(i),F(i)),(1.-S(i))),Z(i));}
  
  // Biomass
  for (i=1;i<=nyrs;i++){
   pred_catch(i) = catage(i) * wt;}

//===================================================================================================
FUNCTION Get_Age_Comp
//===================================================================================================
  
  // Survey
  for (i=1;i<=nyrs;i++){
   for (j=1;j<=nages;j++){
    eac_srv(i,j) = n_srv(i,j)/sum(n_srv(i));}}

  // Fishery
  for (i=1;i<=nyrs;i++){
   for (j=1;j<=nages;j++){
    eac_fish(i,j) = catage(i,j)/sum(catage(i));}}

//===================================================================================================
FUNCTION Evaluate_Objective_Function 
//===================================================================================================

//Objective funtion only for calibration 
  //obj_fun += square(1-1);

// Random effects prior ~N(0,1)
  obj_fun = 0.5*norm2(M_devs);

// Calculate likelihood for survey biomass
  srv_like =  1/(2*square(obs_srv_biom_CV))*sum(square(log(obs_srv_biom)-log(pred_srv)));
         
// Calculate likelihood for catch biomass
  catch_like = 1/(2*square(obs_fish_biom_CV))*sum(square(log(obs_fish_biom)-log(pred_catch)));

// Calculate likelihood for survey age comp
  srv_age_like = -sum(elem_prod(nsamples_srv_age * (obs_ac_srv),log(eac_srv)));

// Calculate likelihood for fishery age comp
  fish_age_like = -sum(elem_prod(nsamples_fish_age * (obs_ac_fish),log(eac_fish)));

// Calculate total likelihood
  obj_fun  += srv_age_like;
  obj_fun  += fish_age_like;
  obj_fun  += srv_like;
  obj_fun  += catch_like;
  
  
//===================================================================================================
FUNCTION write_base_results
//===================================================================================================

  ofstream results("iteration_base.rep");

// With M_devs
  results<<"Convergence TB_end SB_end mean_rec log_M_0 Beta sigma_M M_devs_1985 M_devs_1986  M_devs_1987  M_devs_1988  M_devs_1989  M_devs_1990  M_devs_1991  M_devs_1992  M_devs_1993  M_devs_1994  M_devs_1995  M_devs_1996  M_devs_1997  M_devs_1998  M_devs_1999  M_devs_2000  M_devs_2001  M_devs_2002  M_devs_2003  M_devs_2004  M_devs_2005  M_devs_2006  M_devs_2007  M_devs_2008  M_devs_2009  M_devs_2010  M_devs_2011  M_devs_2012  M_devs_2013  M_devs_2014 M_1985 M_1986  M_1987  M_1988  M_1989  M_1990  M_1991  M_1992  M_1993  M_1994  M_1995  M_1996  M_1997  M_1998  M_1999  M_2000  M_2001  M_2002  M_2003  M_2004  M_2005  M_2006  M_2007  M_2008  M_2009  M_2010  M_2011  M_2012  M_2013  M_2014 mean_F fish_sel_age_2 fish_sel_age_3 fish_sel_age_4 fish_sel_age_5 fish_sel_age_6 fish_sel_age_7 fish_sel_age_8 fish_sel_age_9 fish_sel_age_10 fish_sel_age_11 fish_sel_age_12 fish_sel_age_13 fish_sel_age_14 fish_sel_age_15 fish_sel_age_16 fish_sel_age_17 fish_sel_age_18 fish_sel_age_19 fish_sel_age_20 fish_sel_age_21 fish_sel_age_22 fish_sel_age_23 fish_sel_age_24 fish_sel_age_25 fish_sel_age_26 fish_sel_age_27 fish_sel_age_28 fish_sel_age_29 fish_sel_age_30 fish_sel_age_31 srv_sel_age_2 srv_sel_age_3 srv_sel_age_4 srv_sel_age_5 srv_sel_age_6 srv_sel_age_7 srv_sel_age_8 srv_sel_age_9 srv_sel_age_10 srv_sel_age_11 srv_sel_age_12 srv_sel_age_13 srv_sel_age_14 srv_sel_age_15 srv_sel_age_16 srv_sel_age_17 srv_sel_age_18 srv_sel_age_19 srv_sel_age_20 srv_sel_age_21 srv_sel_age_22 srv_sel_age_23 srv_sel_age_24 srv_sel_age_25 srv_sel_age_26 srv_sel_age_27 srv_sel_age_28 srv_sel_age_29 srv_sel_age_30 srv_sel_age_31 logR rec_devs_1985 rec_devs_1986 rec_devs_1987 rec_devs_1988 rec_devs_1989 rec_devs_1990 rec_devs_1991 rec_devs_1992 rec_devs_1993 rec_devs_1994 rec_devs_1995 rec_devs_1996 rec_devs_1997 rec_devs_1998 rec_devs_1999 rec_devs_2000 rec_devs_2001 rec_devs_2002 rec_devs_2003 rec_devs_2004 rec_devs_2005 rec_devs_2006 rec_devs_2007 rec_devs_2008 rec_devs_2009 rec_devs_2010 rec_devs_2011 rec_devs_2012 rec_devs_2013 rec_devs_2014 init_devs_3 init_devs_4 init_devs_5 init_devs_6 init_devs_7 init_devs_8 init_devs_9 init_devs_10 init_devs_11 init_devs_12 init_devs_13 init_devs_14 init_devs_15 init_devs_16 init_devs_17 init_devs_18 init_devs_19 init_devs_20 init_devs_21 init_devs_22 init_devs_23 init_devs_24 init_devs_25 init_devs_26 init_devs_27 init_devs_28 init_devs_29 init_devs_30 init_devs_31 log_avg_F F_devs_1985 F_devs_1986 F_devs_1987 F_devs_1988 F_devs_1989 F_devs_1990 F_devs_1991 F_devs_1992 F_devs_1993 F_devs_1994 F_devs_1995 F_devs_1996 F_devs_1997 F_devs_1998 F_devs_1999 F_devs_2000 F_devs_2001 F_devs_2002 F_devs_2003 F_devs_2004 F_devs_2005 F_devs_2006 F_devs_2007 F_devs_2008 F_devs_2009 F_devs_2010 F_devs_2011 F_devs_2012 F_devs_2013 F_devs_2014 Fmort_1985 Fmort_1986 Fmort_1987 Fmort_1988 Fmort_1989 Fmort_1990 Fmort_1991 Fmort_1992 Fmort_1993 Fmort_1994 Fmort_1995 Fmort_1996 Fmort_1997 Fmort_1998 Fmort_1999 Fmort_2000 Fmort_2001 Fmort_2002 Fmort_2003 Fmort_2004 Fmort_2005 Fmort_2006 Fmort_2007 Fmort_2008 Fmort_2009 Fmort_2010 Fmort_2011 Fmort_2012 Fmort_2013 Fmort_2014 log_q_srv a50_srv  delta_srv  a50_fish delta_fish tot_biom_1985 tot_biom_1986 tot_biom_1987 tot_biom_1988 tot_biom_1989 tot_biom_1990 tot_biom_1991 tot_biom_1992 tot_biom_1993 tot_biom_1994 tot_biom_1995 tot_biom_1996 tot_biom_1997 tot_biom_1998 tot_biom_1999 tot_biom_2000 tot_biom_2001 tot_biom_2002 tot_biom_2003 tot_biom_2004 tot_biom_2005 tot_biom_2006 tot_biom_2007 tot_biom_2008 tot_biom_2009 tot_biom_2010 tot_biom_2011 tot_biom_2012 tot_biom_2013 tot_biom_2014 spawn_biom_1985 spawn_biom_1986 spawn_biom_1987 spawn_biom_1988 spawn_biom_1989 spawn_biom_1990 spawn_biom_1991 spawn_biom_1992 spawn_biom_1993 spawn_biom_1994 spawn_biom_1995 spawn_biom_1996 spawn_biom_1997 spawn_biom_1998 spawn_biom_1999 spawn_biom_2000 spawn_biom_2001 spawn_biom_2002 spawn_biom_2003 spawn_biom_2004 spawn_biom_2005 spawn_biom_2006 spawn_biom_2007 spawn_biom_2008 spawn_biom_2009 spawn_biom_2010 spawn_biom_2011 spawn_biom_2012 spawn_biom_2013 spawn_biom_2014 nyrs nages spawn_fract recage wt_mat_2  wt_mat_3  wt_mat_4  wt_mat_5  wt_mat_6  wt_mat_7  wt_mat_8  wt_mat_9  wt_mat_10 wt_mat_11 wt_mat_12 wt_mat_13 wt_mat_14 wt_mat_15 wt_mat_16 wt_mat_17 wt_mat_18 wt_mat_19 wt_mat_20 wt_mat_21 wt_mat_22 wt_mat_23 wt_mat_24 wt_mat_25 wt_mat_26 wt_mat_27 wt_mat_28 wt_mat_29 wt_mat_30 wt_mat_31 obs_srv_biom_CV obs_srv_biom_1985 obs_srv_biom_1986 obs_srv_biom_1987 obs_srv_biom_1988 obs_srv_biom_1989 obs_srv_biom_1990 obs_srv_biom_1991 obs_srv_biom_1992 obs_srv_biom_1993 obs_srv_biom_1994 obs_srv_biom_1995 obs_srv_biom_1996 obs_srv_biom_1997 obs_srv_biom_1998 obs_srv_biom_1999 obs_srv_biom_2000 obs_srv_biom_2001 obs_srv_biom_2002 obs_srv_biom_2003 obs_srv_biom_2004 obs_srv_biom_2005 obs_srv_biom_2006 obs_srv_biom_2007 obs_srv_biom_2008 obs_srv_biom_2009 obs_srv_biom_2010 obs_srv_biom_2011 obs_srv_biom_2012 obs_srv_biom_2013 obs_srv_biom_2014 pred_srv_1985 pred_srv_1986 pred_srv_1987 pred_srv_1988 pred_srv_1989 pred_srv_1990 pred_srv_1991 pred_srv_1992 pred_srv_1993 pred_srv_1994 pred_srv_1995 pred_srv_1996 pred_srv_1997 pred_srv_1998 pred_srv_1999 pred_srv_2000 pred_srv_2001 pred_srv_2002 pred_srv_2003 pred_srv_2004 pred_srv_2005 pred_srv_2006 pred_srv_2007 pred_srv_2008 pred_srv_2009 pred_srv_2010 pred_srv_2011 pred_srv_2012 pred_srv_2013 pred_srv_2014 obs_fish_biom_CV obs_fish_biom_1985 obs_fish_biom_1986  obs_fish_biom_1987  obs_fish_biom_1988  obs_fish_biom_1989  obs_fish_biom_1990  obs_fish_biom_1991  obs_fish_biom_1992  obs_fish_biom_1993  obs_fish_biom_1994  obs_fish_biom_1995  obs_fish_biom_1996  obs_fish_biom_1997  obs_fish_biom_1998  obs_fish_biom_1999  obs_fish_biom_2000  obs_fish_biom_2001  obs_fish_biom_2002  obs_fish_biom_2003  obs_fish_biom_2004  obs_fish_biom_2005  obs_fish_biom_2006  obs_fish_biom_2007  obs_fish_biom_2008  obs_fish_biom_2009  obs_fish_biom_2010  obs_fish_biom_2011  obs_fish_biom_2012  obs_fish_biom_2013  obs_fish_biom_2014 pred_catch_1985 pred_catch_1986 pred_catch_1987 pred_catch_1988 pred_catch_1989 pred_catch_1990 pred_catch_1991 pred_catch_1992 pred_catch_1993 pred_catch_1994 pred_catch_1995 pred_catch_1996 pred_catch_1997 pred_catch_1998 pred_catch_1999 pred_catch_2000 pred_catch_2001 pred_catch_2002 pred_catch_2003 pred_catch_2004 pred_catch_2005 pred_catch_2006 pred_catch_2007 pred_catch_2008 pred_catch_2009 pred_catch_2010 pred_catch_2011 pred_catch_2012 pred_catch_2013 pred_catch_2014 nsamples_srv_age nsamples_fish_age srv_like catch_like srv_age_like fish_age_like obj_fun maxgrad"<<endl;
  results<<42<<" "<<tot_biom(nyrs)<<" "<<spawn_biom(nyrs)<<" "<<exp(logR)<<" "<<log_M_0<<" "<<Beta<<" "<<sigma_M<<" "<<M_devs<<" "<<M<<" "<<exp(log_avg_F)<<" "<<fish_sel<<" "<<srv_sel<<" "<<logR<<" "<<rec_devs<<" "<<init_devs<<" "<<log_avg_F<<" "<<F_devs<<" "<<Fmort<<" "<<log_q_srv<<" "<<a50_srv<<" "<<delta_srv<<" "<<a50_fish<<" "<<delta_fish<<" "<<tot_biom<<" "<<spawn_biom<<" "<<nyrs<<" "<<nages<<" "<<spawn_fract<<" "<<recage<<" "<<wt_mat<<" "<<obs_srv_biom_CV<<" "<<obs_srv_biom<<" "<<pred_srv<<" "<<obs_fish_biom_CV<<" "<<obs_fish_biom<<" "<<pred_catch<<" "<<nsamples_srv_age<<" "<<nsamples_fish_age<<" "<<srv_like<<" "<<catch_like<<" "<<srv_age_like<<" "<<fish_age_like<<" "<<obj_fun<<" "<<objective_function_value::pobjfun->gmax<<" "<<endl;

//===================================================================================================
FUNCTION write_mcmc_results
//===================================================================================================
  if(header==1){
   mcmc_results << "M_0,obj_fun" << endl;
   header=0;}
   mcmc_results << M_0 << "," << obj_fun << endl;

//===================================================================================================
REPORT_SECTION
//===================================================================================================
  report<<"obj_fun"<<endl;
  report<<obj_fun<<endl; 
  report<<""<<endl;

  report<<"obs_srv_biom"<<endl;
  report<<obs_srv_biom<<endl; 
  report<<"pred_srv"<<endl;
  report<<pred_srv<<endl; 
  report<<""<<endl;

  report<<"obs_fish_biom"<<endl;
  report<<obs_fish_biom<<endl; 
  report<<"pred_catch"<<endl;
  report<<pred_catch<<endl; 
  report<<""<<endl; 

  report<<"obs_ac_srv"<<endl;
  report<<obs_ac_srv<<endl; 
  report<<"eac_srv"<<endl;
  report<<eac_srv<<endl; 
  report<<""<<endl;

  report<<"obs_ac_fish"<<endl;
  report<<obs_ac_fish<<endl; 
  report<<"eac_fish"<<endl;
  report<<eac_fish<<endl; 
  report<<""<<endl;

  report<<"natage"<<endl;
  report<<natage<<endl; 
  report<<"Z"<<endl;
  report<<Z<<endl; 
