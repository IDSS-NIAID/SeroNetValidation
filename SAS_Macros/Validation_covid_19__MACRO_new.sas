/*********************************************************************************************************************/
/*********************************************validation for covid19 spike_ and nucleocapsid***************************/
/*********************************************************************************************************************/
/***	Before import dataset I should check the variable names, especially the main outcome titer variable, *********/
/***   	I will keep it's name as "Acon", and remove all non_numeric values from it 							 *********/
/***	For LLOQ and ULOQ data   the variable "result" means the calculate result from standardcurve		 *********/
/*********************************************************************************************************************/



%macro importdat(datfile, dat, sheet);
proc import datafile= &datfile
out = &dat 
dbms = xlsx
replace;
sheet = "&sheet";
run;
/*
data  &dat ;
set  &dat ;
*acon =input(acon, best32.);
if assay ="" then delete;
if acon=. then delete;
lg_acon = log(acon);
run;
*/
%mend;

/*

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\COVID19_VAL_Nucleocapsid_Summary070721.xlsx", sheet=LLOQ);

*/



*  %let dat=ULOQ;
%macro lloq(dat);   /* LLOQ and ULOQ will be calculated with this MACRO */
data &dat ;
set &dat;
/*
if var13 ="Range?" then var13 ="";
acon = input(var13, best32.);
acon_st = input(var12, best32.);
*/
sample_idd = substr(sample_id, 1,5);
lg_acon = log(acon);

run;


proc sql noprint;
create table geo_summary_all as
select unique assay, Sample_idd as sample_id,  Dil_Factor,count(*) as num_repplicate, exp(mean(lg_acon)) as geomean,  100*std(acon)/mean(acon) as cv, 
				exp(mean(log(result))) as mean_result, exp(mean(lg_acon))/ Dil_Factor as geomean_conc /* analyst, day ,*/
from &dat
group by  assay, Sample_idd,  Dil_Factor
;
quit;



/* try use the middle three concentration to calculate the theoretical mean */

data &dat._middle_dil;
set &dat;
where 0.5 <= OD <=2.2;
run;
 
proc sql noprint;
create table geo_summary_theo as
select unique assay,  sample_idd as sample_id,   exp(mean(lg_acon)) as theo_mean
from  &dat._middle_dil
group by  assay, Sample_idd 
;
quit;


data geo_summary_all;
merge geo_summary_all geo_summary_theo;
by assay sample_id;
pct_err= 100*abs(theo_mean- geomean)/theo_mean;
run;


/* find the last row which meet the criateria; */
data &dat._final0;
set geo_summary_all;
where .< pct_err<=50 and .< CV <=30;
run;


%if &dat = LLOQ %then %do;
data &dat._final1;
set &dat._final0;
by assay sample_id;
    if last.sample_id ;
	run;
%end;
%if &dat=ULOQ %then %do;
data &dat._final1;
set &dat._final0;
by assay sample_id;
    if first.sample_id ;
	run;
%end;



proc sql noprint;
create table &dat._sum as
select unique assay, mean(geomean_conc) as mean_conc, median(geomean_conc) as median ,std(geomean_conc) as SD ,
		min(geomean_conc) as min, max(geomean_conc) as max, mean(geomean_conc)-1.96*std(geomean_conc) as LCL, mean(geomean_conc)+1.96*std(geomean_conc) as UCL,  
		(mean(geomean_conc)+1.96*std(geomean_conc))*50 as LLOQ, 50*(mean(geomean_conc))  as ULOQ  /* LLOQ is defined as 95% upper limit of the ddiftribution(97.5 percentile) ;  ULOQ will use the mean of ULOQ data*/
from  &dat._final1
group by  assay;
quit;


%mend lloq;



/********************************************************************************************************************************************/
/****************************************************Linearityity data**************************************************************************/
/********************************************************************************************************************************************/

*    %let dat=  linearity; 

  %let  OD_lo1 =0.3;
 %let  OD_up1 =2.0;

  %let  OD_lo2 =0.2;
   %let  OD_up2 =3.0;




%macro linear(dat, OD_lo_1,OD_lo_2, OD_up1, OD_up2);

proc sort data=&dat; by  assay Sample_id analyst day; run;

/**********************************************************************************************/
/***** set range to 0.4<=OD <=2.1***This range was found with mixed model from R code *********/
/***** O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\Chunming_result\Sars_spike\Linearityity_range.R  ****/
/**********************************************************************************************/  
 

data &dat._range; 
set &dat;
if Dil_Factor =50 then delete;
if substr(assay, length(assay)-2 , 3) = "IgG" and (OD > &OD_up1  or OD < &OD_lo1) then delete;
if substr(assay, length(assay)-2 , 3) = "IgM" and (OD > &OD_up2  or OD < &OD_lo2) then delete;
run;


proc sql;
create table geo_sum_&dat._rng as
select unique assay, sample_id,  Dil_Factor, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean,  100*std(acon)/mean(acon) as cv, 
				mean(result) as acon_st_m, exp(mean(lg_acon))/ Dil_Factor as geomean_conc /* analyst, day ,*/
from  &dat._range
group by  assay, Sample_id,  Dil_Factor
;
quit;


proc sql;
create table &dat._range as
select unique assay, sample_id,  exp(mean(log(geomean))) as geomean, min(geomean_conc)*50 as low_limit, max(geomean_conc)*50 as upper_limit
from  geo_sum_&dat._rng 
group by  assay, Sample_id
;
quit;



proc sql;
create table &dat._range_sum as
select unique assay,    median(low_limit) as low_linear_conc, Std(low_limit) as SD_Low, median(upper_limit) as upper_linear_conc ,Std(upper_limit) as SD_up
from  &dat._range 
group by  assay
;
quit;

%mend;


/********************************************************************************************************************************************/
/****************************************************  Cut point   **************************************************************************/
/********************************************************************************************************************************************/

/*  %let dat=cutpoint;  				*/

%macro cut(dat);

proc sql;
create table geo_summary_&dat as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(acon)/mean(acon) as cv 
				 /* analyst, day ,  Sample_id,*/
from &dat
/* where lg_acon^=.  */
group by assay, Sample_id;
quit;


proc sql;
create table outlier as
select unique assay, median(geomean)+3*std(geomean) as upper_lim
from geo_summary_&dat 
group by  assay;
quit;


data combined;   * remove outliers;
merge geo_summary_&dat outlier;
by  assay;
if geomean > upper_lim then delete;
run;


proc sql;
create table &dat._sum as
select unique  assay,count(*) as num_sample, mean(geomean) as mean, median(geomean) as median, min(geomean) as min, max(geomean) as max, std(geomean) as SD,  mean(geomean)+1.96*std(geomean) as &dat
from combined
group by assay;
quit;

%mend;


/* the cutpoint is high, it is better to combine with precision data then do a ROC curve analysis,*/ 


/********************************************************************************************************************************************/
/****************************************************  sensitivity **************************************************************************/
/********************************************************************************************************************************************/

/*         %let dat=Lloq_challenge;   */


%macro sens(dat);

proc sort data=&dat;
by assay sample_id;
run; 


proc sql;
create table geo_summary_sens as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from &dat
/* where lg_acon^=.  */
group by assay, Sample_id;
quit;



data geo_summary_sens1;
set geo_summary_sens;
by assay Sample_id;
assigned_geomean = exp(log(lag(geomean/2)));
if first.assay then assigned_geomean =geomean;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;
run;

data geo_summary_sens2;
set geo_summary_sens1;
by assay Sample_id;
where .< pct_err<=50 and .< CV <=30;
if last.assay ;
run;

%mend;



/********************************************************************************************************************************************/
/****************************************************  accuracy **************************************************************************/
/********************************************************************************************************************************************/


/*
%let dat=ACCURACY;

*/


%macro acc(dat);


proc sort data=&dat;
by assay Analyst sample_id;
run; 

/****************by analyst *************************/
proc sql;
create table geo_analsyt_&dat as
select unique assay, Analyst, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from &dat
/* where lg_acon^=.  */
group by assay,   Analyst, Sample_id;
quit;


data geo_analsyt_acc_IgG;
set geo_analsyt_&dat;
where substr(assay, length(assay)-2, 3)="IgG";
by assay Analyst Sample_id;
assigned_geomean = lag(geomean/3);  * how to set up the true value for each dilution? ; 
if  first.analyst then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;


data geo_analsyt_acc_IgM;
set geo_analsyt_&dat;
where substr(assay, length(assay)-2, 3)="IgM";
by assay Analyst Sample_id;
assigned_geomean = lag(geomean/2);  * how to set up the true value for each dilution? ; 
if  first.analyst then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;


data combined_analyst;
set geo_analsyt_acc_IgG geo_analsyt_acc_IgM;
run;

/************************  overall ************************************************/

proc sql;
create table geo_ovarall_&dat as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from &dat 
/* where lg_acon^=.  */
group by assay,    Sample_id;
quit;


data geo_all_IgG;
set geo_ovarall_&dat;
where substr(assay, length(assay)-2, 3)="IgG";
by assay  Sample_id;
assigned_geomean = lag(geomean/3);  * how to set up the true value for each dilution? ; 
if  first. assay  then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;

data geo_all_IgM;
set geo_ovarall_&dat;
where substr(assay, length(assay)-2, 3)="IgM";
by assay  Sample_id;
assigned_geomean = lag(geomean/2);  * how to set up the true value for each dilution? ; 
if  first. assay  then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;


data comnined_all;
set geo_all_IgG geo_all_IgM;
run;

%mend;

/********************************************************************************************************************************************/
/****************************************************presicion data**************************************************************************/

/*   %let dat=Precision;			*/

%macro pres(dat);

data pres;
set &dat ;
day_c=day;
*assay=HPV_type;
run;


proc sort data=pres;
by assay sample_id analyst day;

run;

proc sql;
create table geo_summary_pres as
select unique  Assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon))/ mean(exp(lg_acon)) as cv
				 /* analyst, day ,*/
from pres
group by Assay, Sample_id;
quit;

proc sort data=pres;
by assay sample_id;
run;

proc sql;
create table  che as
select unique assay, Analyst, day, count(*)
from pres
group by assay, Analyst, day;
quit;



data pres;
set pres;
if assay="" then delete;
run;


%include "O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\mixed_model for_ICC_CV_July2021.sas";  * start to call macro to calculate CV;
 

%ICC_CV(pres); 



%mend;





/********************************************************************************************************************************************/
/****************************************************Carry over**************************************************************************/



%macro  carry(dat);

proc sql;
create table geo_summary_carry as
select unique  Assay, Sample_id, analyst, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,*/
from CARRYOVER
group by Assay, Sample_id;
quit;

%mend; 

/********************************************************************************************************************************************/
/****************************************************Specificity**************************************************************************/

/*    %let dat=Specificity;   */

%macro spec(dat);
data Specificity;
set &dat;
*assay=HPV_type;
antibody= substr(sample_id, (prxmatch("/PC/", Sample_ID )),3);
run;
quit;

/* find the unspiked data  */



proc sort data=Specificity;
by assay Spiked_Antigen Sample_ID;
run;

proc sql;
create table Specificity_sum as 
select unique  Assay, antibody,Spiked_Antigen, Sample_ID, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Specificity
group by Assay, antibody, Spiked_Antigen,Sample_ID;
quit;

data unspiked;
set  Specificity_sum;
where Spiked_Antigen="Unspiked";
assigned=geomean;
keep Assay  antibody  assigned;
run;




data  Specificity_sum1;
merge Specificity_sum unspiked;;
by Assay  antibody ;
pct_err = 100*(geomean-assigned )/assigned;
run;


proc sort data= Specificity_sum1;
by  Assay antibody  descending Spiked_Antigen; 
run:



data summ;
set Specificity_sum1;
*where  pct_err1 ^=0;
format geomean assigned 8.2 CV 8.4 pct_err 8.2;
drop antibody ;
run;
%mend;


/********************************************************************************************************************************************/
/****************************************************stability ******************************************************************************/
/********************************************************************************************************************************************/
/********************************** freeze and thraw*****************************************************************************************/



/*   %let dat=STABILITY_Heat;  */  

%macro stable(dat);

proc sort data=&dat;
by assay sample_id;
run;

%if &dat=STABILITY_FREEZE_THAW %then %do;
	data Stable;
	set &dat;
	where SAmple_id ^="";
	if sample_id in ("Low 1x" "High 1x") then trt=0;
	else trt=1;
	grp=substr(sample_id, 1,3); 
	run;
%end;

%else %if &dat=STABILITY_HEAT %then %do;
 	data Stable;
	set &dat;
		if Treatment="Heat" then Trt=1;
		else trt=0; 

		if substr(sample_id, 1,4)="Heat" then grp=substr(sample_id, 6,3); 
		
	 	else grp=substr(sample_id, 1,3);
		
	run;
%end;

%else %do;

	data Stable;
	set &dat;
	where SAmple_id ^="";
	if Treatment="" then Trt=0; 
	else trt=1;
	grp=substr(sample_id, 1,3); 
	run;
%end;

data Stable1;
set Stable;
 
if grp in ("HIG" "Hig") then grp ="HIG";
if grp in ("MED" "Med") then grp ="MED";
if grp in ("Low" "LOW") then grp ="LOW";
run;
proc sort data=Stable;
by assay sample_id;
run;

proc mixed data=Stable1;
class assay grp trt(ref='0');
by assay;
model lg_acon =  trt;
random int/subject=grp;
run;


proc sql;
create table Stable_trt_sum as 
select unique  Assay, sample_id, grp, trt, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable1
group by Assay, grp, trt;
quit;

proc sort data=Stable_trt_sum;
by assay grp trt;
run;

data Stable_trt_sum;
set Stable_trt_sum;
by  assay grp trt;;
val=lag(geomean);
if first.grp then val=geomean;
pct_err = 100*(geomean-val)/val;
run;
%mend;




/*************************************************************************************************************************************/
/**********************************Lot to Lot*****************************************************************************************/
/*************************************************************************************************************************************/

/*******************************************lot to lot Antigen*****************************************************************************************/
/*        %let dat=lot_a;											*/

%macro lot(dat);
proc sort data= &dat;
by Assay sample_id;
run;

data Lot_A;
set &dat;
 lg_acon=log(acon);
 lot=Lot_Status__Old_or_New;
run;
title;



proc mixed data=Lot_A method=reml;
class Assay sample_id Lot;
by Assay ;
model lg_acon=Assay  Lot/s ddfm=res;
random int / subject= sample_id;
run;

data lot_A_N lot_A_O;
set  Lot_A;
if lot="New" then output lot_A_N;
if lot in ("old" "Old") then output lot_A_O;
run;

proc sql;
create table combine_wide as
select a.*, b.lg_acon as lg_acon_new, b.acon as acon_new
from  lot_a_O a left join lot_a_N b on
a.assay=b.assay and a.sample_id =b.sample_id
;
quit;

data  combine_wide1;
set  combine_wide;
pct_err = 100 * (acon_new - acon )/acon;
run;

proc sql noprint ;
create table nopass as
select unique assay , count(*) as number_no_pass
from Combine_wide1
where abs(pct_err) >=25
group by assay ;
quit;

proc sql noprint ;;
create table summary as
select unique assay,  count(*) as Num_exp, mean(pct_err) as mean_pct_err, Std(pct_err) as SD_pct_err, min(pct_err) as min_pct_err, max(pct_err) as max_pct_err
from Combine_wide1
group by assay ;
quit;

data summ;
merge summary nopass;
by assay;
if number_no_pass=. then number_no_pass=0;
num_pass =Num_exp-number_no_pass;
pass_pct =100*(num_pass/Num_exp);
run;

%mend;

