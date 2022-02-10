/********************************program: ***************************************************************************************************/  
/* O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\ELISA_Validation\Validation_Analysis_covid_19nucleocapsid_try0927.sas      *******/ 
/* O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\ELISA_Validation\Validation_Analysis_covid_19nucleocapsid_try0927.sas      *******/ 
 
/****Data Folder: O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA
/***Data: 	COVID19_VAL_Spike_Summary_vFinal_23SEP21.xlxs
`		    COVID19_VAL_Nucleocapsid_Summary_Final_23SEP21.xlsx

***********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************validation for covid19 spike_*******************************************/
/*********************************************************************************************************************/

/***************************************************************************LLOQ data**********************************************************************/

%include "O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\MACROs\Validation_covid_19__MACRO_new.sas";

%let file = COVID19_VAL_Spike_Summary_vFinal_23SEP21;

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=LLOQ, sheet=LLOQ);

data lloq;
set lloq;
sample_idd = substr(sample_id, 1,5);
drop sample_id;
run;

/***************************************************************************ULOQ data***********************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=ULOQ, sheet=ULOQ);

data uloq;
set uloq;
sample_idd = substr(sample_id, 1,5);
drop sample_id ;

run;

/****************************************************Linearity data**************************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Linear, sheet=LINEARITY);

data linear;
set linear;
sample_idd = substr(sample_id, 1,5);
drop sample_id ;
run;

proc sort data=linear; by  assay Sample_idd analyst day; run;

/*******************  combine LLOQ ULOQ Linearity  *********************************/

data all_linear;
 set lloq linear Uloq;
 if acon=. then delete;
 run;

proc sql; 
 select min(OD), max(OD) into: ODmin  , :ODmax  
 from  all_linear;
 quit;


data all_linear_1;
 set all_linear; /*  where lg_OD=.; run;  */
 lg_OD = log((OD- &ODmin )/(&ODmax- OD));
lg_dil=log(Dil_Factor);
 run;


/*********************************************************************************************************/
/***** set range to for IgG and IgM separately This range was found with mixed model from R code *********/
/***** O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\Chunming_result\Sars_spike\linearity_range.R  ***************/
/*********************************************************************************************************/  
 
data linear_range;
set all_linear;

if Assay ="Spike IgG" and ( 0.12 > OD or OD > 2.9 ) then delete;
if Assay ="Spike IgM" and ( 0.15 > OD or OD > 2.7 ) then delete;
run;

proc sql;
create table geo_sum_linear_rng as
select unique assay, Sample_idd as sample_id,  Dil_Factor, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean,  100*std(acon)/mean(acon) as cv, 
				mean(result) as acon_st_m, exp(mean(lg_acon))/ Dil_Factor as geomean_conc /* analyst, day ,*/
from linear_range
group by  assay, Sample_idd,  Dil_Factor
;
quit;

proc sql;
create table lin_range_sum as
select unique assay, sample_id,  exp(mean(log(geomean))) as geomean, min(geomean_conc)*50 as low_limit, max(geomean_conc)*50 as upper_limit
from  geo_sum_linear_rng 
group by  assay, Sample_id
;
quit;

data lin_range_sum;
set lin_range_sum;
if low_limit =upper_limit then delete;
run;


proc sql;
create table lin_range_final as
select unique assay, count(*) as num_sample, median(low_limit) as median_lo , mean(low_limit) as mean_lo, std(low_limit)/sqrt(count(*)) as SE_lo, 
									 median(upper_limit) as median_up , mean(upper_limit) as mean_up, std(upper_limit)/sqrt(count(*)) as SE_up 

from lin_range_sum
group by assay;
quit;



proc sort data=linear_range;
by assay sample_idd od;
run;


/***************************************************************************************************************************************************************************************************************************************************************/

proc sql;
create table geo_summary_linear as
select unique assay, Sample_idd as sample_id,  Dil_Factor, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, median(OD) as median_OD, 100*std(acon)/mean(acon) as cv, 
				exp(mean(log(acon_st))) as acon_st_m, exp(mean(lg_acon))/ Dil_Factor as geomean_conc /* analyst, day ,*/
from linear
group by  assay, Sample_idd,  Dil_Factor
;
quit;

data linear_theo1;
set linear;
where Dil_Factor =1350;
run;


proc sql;
create table linear_theo as
select unique assay, Sample_idd as sample_id, exp(mean(lg_acon)) as theo_mean
from linear_theo1
group by  assay, Sample_id;
quit;


data geo_summary_lin_all;
merge geo_summary_linear  linear_theo;;
by assay Sample_id;
pct_err= 100*abs(theo_mean- geomean)/theo_mean;
run;


data linear_CV_pece;   
set geo_summary_lin_all;
where .< pct_err<=50 and .< CV <=30;
run;


proc mixed data=geo_summary_lin_all;
class assay sample_id; 
by Assay;
model median_OD=acon_st_m/s;
random intercept /subject=Sample_id;
run;

proc corr data =geo_summary_lin_all nomiss outp=CorrOutp(type=corr);
by assay sample_id;
var median_OD  geomean_conc;
*ods output PearsonCorr=gemma.;
run;

/*********************************************************************************************/

data linear_range;
set all_linear;
where  0.05 <= OD <= 2.10;
run;

proc sql;
create table assigned as
select unique assay,sample_idd as sample_id, exp(mean(lg_acon)) as assigned_mean
from   linear_range
group by  assay ,sample_idd;
quit;


proc sql;
create table LOQ_sum as
select unique assay,sample_idd as sample_id, Dil_Factor, exp(mean(lg_acon)) as geomena, mean(result) as concentration, std(result)/mean(result)as CV_con, 100*std(acon)/mean(acon)as CV 
		
from   all_linear
group by  assay ,sample_idd, Dil_Factor;
quit;


data LOQ_sum2;
merge LOQ_sum assigned ;
by assay sample_id;
pct_err = abs(100*(geomena - assigned_mean))/assigned_mean;
run;


/*********************************************************************************************************/
/************   use combined sample to find lloq and ULOQ  ***********************************************/
/************ for lloq, only use those with low limits(a.k.a. there is large CV/pct_err in lower end)**/
/************ for ULOQ, only use those with upper limits(a.k.a. there is large CV/pct_err in hi end)  **/
/*********************************************************************************************************/

/* find sample for LLOQ and ULOQ; */
proc sql;
create table lloq_id as
select unique assay, sample_id
from LOQ_sum2
/* where Dil_Factor>= 4050  and (pct_err > 50 or CV >30) */
group by assay, sample_id;
quit;


proc sql;
 create table LLOQ_sum as
 select a.*, b.*
 from lloq_id a left join LOQ_sum2 b on 
 a.assay=b.assay and a.sample_id=b.sample_id;
 quit;



proc sort data=LLOQ_sum;
by assay sample_id Dil_Factor; run;

data lloq_sum2;
set LLOQ_sum;
by assay sample_id Dil_Factor;
where pct_err <= 50 and CV <=30;
if last.sample_id;
*if Dil_Factor<=1350 then delete;
run;

/*************check LgM LIN_2 and LIN_4 ******************/

data Lin_2_4;
set Linear;
where assay="Spike IgM" and sample_idd in ("LIN_2" "LIN_4");
if acon=. then delete;
run;


proc sort data=Lin_2_4; by sample_idd day Dil_Factor;
run;

/******************************************************************/
/******* remove outlier  Q3 +_ 2IQR ***************************/

proc means data=lloq_sum2 n median std p25 p75 p95;
class assay;
var  concentration;
run;

/* find the outliers bound */
/*
Assay 		N Obs 	N   Median      Std_Dev 25th Pctl 75th Pctl 95th Pctl 
Spike IgG 	20 		20 0.3394107 0.2301043  0.2093929 0.4026310 0.8380000 
Spike IgM 	20 		20 0.1369095 4.0841413  0.0792250 0.2282679 13.4031667 


IQR_IgG = 0.4026310  -0.2093929 = 0.1932381;
IQR_IgM = 0.2282679 -0.0792250  = 0.1490429;

uplimit_IgG = Q3 +3*IQR = 0.4026310 + 3*0.1932381 =0.9824;
uplimit_IgM = Q3 +3*IQR = 0.2282679 + 3*0.1490429 =  0.6754

*/

/**************************************************/

data lloq_sum3 ;
set lloq_sum2 ;
if assay = "Spike IgG" and concentration > 0.9824 then delete;
if assay = "Spike IgM" and concentration > 0.6754 then delete;
run;

proc sql;
create table lloq_final as
select unique assay, count(*) as N,  mean(concentration) as mean, median(concentration) as median, std(concentration) as SD, 50*(mean(concentration)+ 1.645*std(concentration)/sqrt(count(*))) as lloq 
from lloq_sum3
group by assay;
quit;


/*******************************************************************************************************************/

proc sql;
create table Uloq_id as
select unique assay, sample_id
from LOQ_sum2
/* where Dil_Factor<= 12150  and (pct_err > 50 or CV >30) */
group by assay, sample_id;
quit;


proc sql;
 create table ULOQ_sum as
 select a.*, b.*
 from uloq_id a left join LOQ_sum2 b on 
 a.assay=b.assay and a.sample_id=b.sample_id;
 quit;


proc sort data=uLOQ_sum;
by assay sample_id Dil_Factor; run;
 
data uloq_sum2;
set uLOQ_sum;
by assay sample_id Dil_Factor;
where pct_err <= 50 and CV <=30;
if first.sample_id;
*if Dil_Factor>=4050 or concentration < 4 then delete;
run;


/******************************************************************/
/******* remove outlier  Medain +_ 2IQR ***************************/

proc means data=uloq_sum2 n median std p25 p75 p95;
class assay;
var  concentration;
run;

/* find the outliers bound */
/*
Assay      N Obs N  Median     Std Dev   25th Pctl   75th Pctl   95th Pctl 
Spike IgG  19 	19  0.3361786 7.9363026 0.1722143   0.5128438    34.9165714 
Spike IgM 20 	20  0.1418023 4.0722689 0.0792250   0.3311607    13.4031667 

IQR_IgG = 26.1102174 -15.0312903= 11.07;
IQR_IgM = 12.7933929 -2.1495000= 10.64;

uplimit_IgG = median - 2*IQR = 22.8680000 - 2*11.07 =1.017;
uplimit_IgM = median - 2*IQR = 0.1418023 + 2*0.2519 =  0.6456

*/

/**************************************************/
proc sql;
create table uloq_final as
select unique assay,count(*) as N,  mean(concentration) as mean, median(concentration) as median, std(concentration) as SD, 50*mean(concentration) as uloq 
from uloq_sum2
group by assay;
quit;


proc univariate data =lloq_uloq_sum;
by assay;
var lloq uloq;
run;
/* extreme values*/
/* Spike IgG (ULOQ2) lloq=34.92; 
   Spike IgG (ULOQ3) uloq= 4.36; Spike IgG (ULOQ4) uloq= 2.61; Spike IgG (LLOQ5) uloq= 9.65;
   Spike IgM (LIN_2) lloq=12.90; Spike IgM (Lin4) lloq=13.91;
   Spike IgM (Lin3)  uloq=0.725; Spike IgM (LLOQ4) lloq=2.15; Spike IgM (LLOQ7) uloq=1.28; Spike IgM (LLOQ1) uloq=2.23; Spike IgM (ULOQ8) uloq=1.62;
*/

data lloq_uloq_sum1;  /* remove outlieters*/
set lloq_uloq_sum;
if uloq <=3 then uloq=.;
if (Assay ="Spike IgG" and sample_id = 'ULOQ2') or  (Assay ="Spike IgM" and sample_id in ('LIN_4' 'LIN_2')) then lloq=.;
if (Assay ="Spike IgG" and sample_id in ('ULOQ3' 'ULOQ4' 'LLOQ5')) or  (Assay ="Spike IgM" and sample_id in ('LIN_1' 'LIN_3' 'LLOQ3' 'LLOQ4' 'LLOQ5' 'LLOQ7' 'LLOQ1' 'ULOQ7' 'ULOQ8')) then uloq=.;
run;

proc means data=lloq_uloq_sum1 n mean median std p95 p99 maxdec=3;
class assay;
var lloq uloq;
run;
 



/********************************************************************************************************************************************/
/*******************************************  Cut point   no updates so no need to re-do*9-28-2021******************************************/
/********************************************************************************************************************************************/


%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=cut, sheet=CUTPOINT);


proc sql;
create table geo_summary_cut as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, exp(mean(lg_acon) 
				 /* analyst, day ,  Sample_id,*/
from cut
/* where lg_acon^=.  */
group by assay, Sample_id;
quit;



proc sql;
create table cutpoint as
select unique  assay,count(*) as num_sample, mean(geomean) as mean, median(geomean) as median, min(geomean) as min, max(geomean) as max, std(geomean) as SD,  mean(geomean)+2*std(geomean) as Cutpoint
from geo_summary_cut 
group by assay;
quit;

/* the cutpoint is high, it is better to use ROC curve analysis*/ 


/********************************************************************************************************************************************/
/******************************************  sensitivity-- no updates so no need to re-do*9-28-2021******************************************/
/********************************************************************************************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=sens, sheet=LLOQ_Challenge);


proc sort data=sens;
by assay sample_id;
run; 


proc sql;
create table geo_summary_sens as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from sens
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



/********************************************************************************************************************************************/
/*******************************************  accuracy -- no updates so no need to re-do*9-28-2021******************************************/
/********************************************************************************************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=acc, sheet=ACCURACY);

proc sort data=acc;
by assay Analyst sample_id;
run; 


proc sql;
create table geo_summary_acc as
select unique assay, Analyst, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from acc
/* where lg_acon^=.  */
group by assay,   Analyst, Sample_id;
quit;


data geo_summary_acc1_IgG;
set geo_summary_acc;
where assay="Spike IgG";
by assay Analyst Sample_id;
assigned_geomean = lag(geomean/3);  * how to set up the true value for each dilution? ; 
if  first.analyst then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;


data geo_summary_acc1_IgM;
set geo_summary_acc;
where assay="Spike IgM";
by assay Analyst Sample_id;
assigned_geomean = lag(geomean/2);  * how to set up the true value for each dilution? ; 
if  first.analyst then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;


data combined;
set geo_summary_acc1_IgG geo_summary_acc1_IgM;
run;


proc sql;
create table geo_summary_acc as
select unique assay, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,  Sample_id,*/
from acc
/* where lg_acon^=.  */
group by assay,    Sample_id;
quit;

data geo_summary_acc1_IgM;
set geo_summary_acc;
where assay="Spike IgM";
by assay  Sample_id;
assigned_geomean = lag(geomean/2);  * how to set up the true value for each dilution? ; 
if  first. assay  then assigned_geomean =geomean;   *the first oone as true value?????   ;
pct_err = 100*abs(geomean-assigned_geomean)/assigned_geomean;

run;

/********************************************************************************************************************************************/
/****************************************presicion data*-- no updates so no need to re-do*9-28-2021******************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=pres, sheet=PRECISION);

data empty;set pres;
where assay="";
run;


data pres;
set pres;
if assay="" then delete;
run;

%include "O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\mixed_model for_ICC_CV_July2021.sas";  * start to call macro to calculate CV;
 

%ICC_CV(pres); 


ODS rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\Chunming_result\CVs and ICC for SARs spike_precision.rtf";
     proc print data=icc noobs split='*';
	 var assay cvwd cvbd cvbt cv icc;
     label cvwd='Within Day CV*' cvbd='Between Day CV*' cvbt='Between Tech CV*' cv='Overall CV*' icc='ICC*';
     title "CVs and ICC for Overall reproducibility test";
	 run;
ODS rtf close;


/********************************************************************************************************************************************/
/********************************************Carry over*-- no updates so no need to re-do*9-28-2021******************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=CARRYOVER, sheet=CARRYOVER);


proc sql;
create table geo_summary_carry as
select unique  Assay, Sample_id, analyst, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
				 /* analyst, day ,*/
from CARRYOVER
group by Assay, Sample_id;
quit;



/********************************************************************************************************************************************/
/****************************************************Specificity**************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Specificity, sheet=Specificity);

data Specificity;
set Specificity;
antibody= substr(sample_id, (prxmatch("/PC/", Sample_ID )),3);
run;


proc sort data==Specificity;
by assay ANTIBODY  sample_id;
run;

proc mixed data=Specificity;
by assay;
class antibody(ref='PC1') Assay Spiked_Antigen;
model lg_acon=assay antibody/s ddfm=res;
random int / subject=Spiked_Antigen;
run;   /* IgG p <0.001 ;  IgM =0.0036;  overall p < 0.0001*/


proc sort data=Specificity;
by assay DESCENDING Spiked_Antigen  Sample_ID;
run;

proc sql;
create table Specificity_sum as 
select unique  Assay, ANTIBODY, Spiked_Antigen, Sample_ID, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Specificity
group by Assay, ANTIBODY, Spiked_Antigen, Sample_ID;
quit;

proc sort data=Specificity_sum;
by assay ANTIBODY DESCENDING Spiked_Antigen ;
run;

data  Specificity_sum1;
set  Specificity_sum;
by Assay ANTIBODY DESCENDING Spiked_Antigen ;
val=lag(geomean);
if first.ANTIBODY then val=geomean;
pct_err = 100*(geomean-val)/val;

run;


/********************************************************************************************************************************************/
/****************************************************stability ******************************************************************************/
/********************************************************************************************************************************************/
/********************************** freeze and thraw*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Stable_frez_thr, sheet=STABILITY_FREEZE-THAW);


data Stable_1;
set Stable_frez_thr;
run;


proc sort data=Stable_1;by assay sample_id;
run;

data Stable_1;
set Stable_1;
where SAmple_id ^="";
if Treatment = "1X" then freeze=1;
else if Treatment = "5X" then freeze=5;
else freeze=10;
run;


proc freq data=Stable_1;
table SAmple_id *Treatment/norow nocol nopercent nocum;
run;

proc mixed data=Stable_1;
class assay  SAmple_id freeze(ref='1');
by assay;
model lg_acon =  freeze;
random int/subject= SAmple_id;
contrast 'Freeze 10x vs 5x' 	freeze 1 -1  0;
contrast 'Freeze 10x vs normal' freeze 1  0 -1;
contrast 'Freeze 5x  vs normal' freeze 0  1 -1;
run;


proc sql;
create table Stable_frz_sum as 
select unique  Assay, sample_id, Treatment, freeze, Treatment, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable_1
group by Assay, sample_id, freeze;
quit;

proc sort data=Stable_frz_sum;
by assay sample_id freeze;
run;

data Stable_frz_sum_first;
set Stable_frz_sum;
by  assay sample_id freeze;
if first.sample_id ;
val=geomean;
keep assay sample_id val;
run;


data Stable_frz_summary;
merge Stable_frz_sum Stable_frz_sum_first;
by assay sample_id;
pct_err = 100*(geomean -val)/val;
run;


/********************************** stability _ Lipedia*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Stable_lipe, sheet=STABILITY_LIPEMIA);

data Stable_2;
set  Stable_lipe;
lg_acon=log(acon);
run;



proc freq data=Stable_2;
table assay*sample_id*Treatment/norow nocol nopercent nocum;
run;

proc mixed data=Stable_2;
class assay Sample_ID Treatment(ref='Normal');
by assay;
model lg_acon =  Treatment/s;
random int/subject=Sample_ID;
run;



proc sql;
create table Stable_lipe_sum as 
select unique  Assay, sample_id, Sample_ID, Treatment, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable_2
group by Assay, Sample_ID, Treatment;
quit;

proc sort data=Stable_lipe_sum ;
by assay Sample_ID descending treatment;
run;

data Stable_lipe_sum ;
set Stable_lipe_sum ;
by  assay Sample_ID descending treatment;
val=lag(geomean);
if first.Sample_ID then val=geomean;
pct_err = 100*(geomean-val)/val;
run;




/**********************************STABILITY_BILIRUBIN*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=STABILITY_BILIRUBIN, sheet=STABILITY_BILIRUBIN);


data Stable_3;
set  STABILITY_BILIRUBIN;
if Treatment="" then Treatment="Normal";
run;



proc freq data=Stable_3;
table assay*sample_id*Treatment/norow nocol nopercent nocum;
run;

proc mixed data=Stable_3 method=REML;
class assay sample_id Treatment(ref="Normal");
by assay;
model lg_acon =  Treatment/s;
random int/subject=sample_id;
contrast 'Bilirubin high vs low' 	Treatment 1 -1  0;
contrast 'Bilirubin high vs normal' Treatment 1  0 -1;
contrast 'Bilirubin low vs normal' 	Treatment 0  1 -1;
run;


proc sql;
create table Stable_BILI_sum as 
select unique  Assay, sample_id, Treatment, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable_3
group by Assay,  sample_id, Treatment;
quit;

proc sort data=Stable_BILI_sum ;
by assay  sample_id descending treatment;
run;

data Stable_BILI_sum ;
set Stable_BILI_sum ;
by  assay sample_id descending treatment;
val=lag(geomean);
if first.sample_id then val=geomean;
pct_err = 100*(geomean-val)/val;
run;



/**********************************STABILITY_HEMOLYSIS*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=STABILITY_HEMOLYSIS, sheet=STABILITY_HEMOLYSIS);

data Stable_4;
set  STABILITY_HEMOLYSIS;
if Treatment="" then Treatment="Normal";
run;


proc freq data=Stable_4;
table assay*sample_id*Treatment/norow nocol nopercent nocum;
run;

proc mixed data=Stable_4;
class assay sample_id Treatment(ref='Normal');
by assay;
model lg_acon =  Treatment/s;
random int/subject=sample_id;
run;


proc sql;
create table Stable_hemo_sum as 
select unique  Assay, sample_id,  Treatment, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable_4
group by Assay, sample_id, Treatment;
quit;

proc sort data=Stable_hemo_sum ;
by assay sample_id descending treatment;
run;

data Stable_hemo_sum ;
set Stable_hemo_sum ;
by  assay sample_id descending treatment;
val=lag(geomean);
if first.sample_id then val=geomean;
pct_err = 100*(geomean-val)/val;
run;


/**********************************STABILITY_HEAT*****************************************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=STABILITY_HEAT, sheet=STABILITY_HEAT);

data Stable_5;
set STABILITY_HEAT;
 if substr(Sample_id,1,3)="Hea" then grp= substr(sample_id,6);
 else grp= substr(sample_id,1);
 lg_acon=log(acon);

run;


proc sort data= Stable_5;by assay grp descending treatment;
run;
proc freq data=Stable_5;
table assay*grp*Treatment/norow nocol nopercent nocum;
run;

proc mixed data=Stable_5;
class assay grp Treatment(ref='Normal');
by assay;
model lg_acon =  Treatment/s;
random int/subject=grp;
run;



proc sql;
create table Stable_heat_sum as 
select unique  Assay, sample_id, grp, Treatment, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean, 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV
from Stable_5
group by Assay, grp, Treatment;
quit;

proc sort data=Stable_heat_sum ;
by assay grp descending treatment;
run;

data Stable_heat_sum ;
set Stable_heat_sum ;
by  assay grp descending treatment;
val=lag(geomean);
if first.grp then val=geomean;
pct_err = 100*(geomean-val)/val;
run;


/*************************************************************************************************************************************/
/**********************************Lot to Lot*****************************************************************************************/
/*************************************************************************************************************************************/

/*******************************************lot to lot Antigen*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Lot_Antigen, sheet=Lot-to-Lot Antigen);


proc sort data= Lot_Antigen;
by Assay sample_id;
run;

data Lot_Antigen;
set Lot_Antigen;
 lg_acon=log(acon);
 lot=Lot_Status__Old_or_New;
run;

proc mixed data=Lot_Antigen method=reml;
class Assay sample_id Lot(ref='Old');
by Assay ;
model lg_acon=Assay  Lot/s ddfm=res;
random int / subject= sample_id;
run;

data lot_antigen_N lot_antigen_O;
set  Lot_Antigen;
if lot="New" then output lot_antigen_N;
if lot="Old" then output lot_antigen_O;
run;

proc sql;
create table combine_wide as
select a.*, b.lg_acon as lg_acon_new, b.acon as acon_new
from  lot_antigen_O a left join lot_antigen_N b on
a.assay=b.assay and a.sample_id =b.sample_id
;
quit;

data  combine_wide1;
set  combine_wide;
pct_err = 100 * (acon_new - acon )/acon;
pct_err_log = 100*(lg_acon_new -lg_acon)/lg_acon;
run;

proc means data=combine_wide1 n mean std stderr min median max maxdec=2;
class assay;
var pct_err;
run;

ODS RTF file="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\Chunming_result\result_28Sep21\COVID19_Spike_lot_antigen.rtf";
proc print data=combine_wide1 noobs split='*';
	 var assay Sample_ID  pct_err;
     label pct_err='percent error between old and new lot*';
 	 run;


ODS rtf Close;

ods html close;

/*******************************************lot to lot Conjugate*****************************************************************************************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\&file", dat=Lot_conjugate, sheet=Lot-to-Lot Conjugate);


proc sort data=  Lot_conjugate;
by Assay sample_id;
run;

data  Lot_conjugate;
set  Lot_conjugate;
if assay="" then delete;
if acon=. then delete;
 lg_acon=log(acon);
 lot=Lot_Status__Old_or_New;
 if lot ="old" then lot ="Old";
run;

proc mixed data=Lot_conjugate method=reml;
class Assay sample_id Lot(ref='Old');
by Assay ;
model lg_acon =Assay Lot/s ddfm=res;
random int / subject= sample_id;
run;


proc sql;
create table summ as
select unique assay,   lot, sample_id,exp(mean(lg_acon)) as geomean, 100*std(acon)/mean(acon) as CV
from Lot_conjugate
group by  assay , lot,  sample_id;
quit;


data lot_conjugate_N lot_conjugate_O;
set  summ;
if lot="New" then output lot_conjugate_N;
if lot="Old" then output lot_conjugate_O;
run;

proc sql;
create table combine_wide_c as
select a.*, b.geomean as geomean_new
from  lot_conjugate_O a left join lot_conjugate_N b on
a.assay=b.assay and a.sample_id =b.sample_id
;
quit;

data  combine_wide_c;
set  combine_wide_c;
pct_err = 100*(geomean_new -geomean)/geomean_new;

run;

proc means data=combine_wide_c n mean std stderr min median max maxdec=2;
class assay;
var  pct_err;
run;



ODS RTF file="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\Chunming_result\result_28Sep21\COVID19_Spike_lot_conjugate.rtf";
proc print data=combine_wide_c noobs split='*';
	 var assay Sample_ID  pct_err;
     label pct_err='percent error between old and new lot*';
 	 run;


ODS rtf Close;


/*****************specificity follow up data ************************************************************/


proc import datafile= 'O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\ELISA\COVID19_VAL_Spike_Summary_vFinal_23SEP21.xlsx'
out = Spec
dbms = xlsx
replace;
sheet = "Specificity Follow-Up";
run;


data spec;
set spec;
if assay="Spike IgG" and acon > 31.22 then positive=1;
else if assay="Spike IgM" and acon > 81.69 then positive=1;
else positive=0;
run;

proc sql;
create table spec_sum as
select unique assay,  count(*)-sum(positive) as true_neg, count(*) as total_num, 100*(1- sum(positive)/count(*)) as sepcificity
from spec
group by assay;
quit;
