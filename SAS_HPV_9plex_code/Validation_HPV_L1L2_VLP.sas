
/*********************************************************************************************************************/
/*********************************************validation for covid19 Spike********************************************/
/*********************************************************************************************************************/
/***	Before import dataset I should check the variable names, especially the main outcome titer variable, *********/
/***   	I will keep it's name as "Acon", and remove all non_numeric values from it 							 *********/
/***	For LLOQ and ULOQ data   the variable "result" means the calculate result from standardcurve		 *********/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/**********MAKE SURE ALL "Range?"  "Masked" "  "(space) are changed to BLANK(but not space)		**********************/
/*********************************************************************************************************************/



/****************************************************/
/* 	 		root directory				*/

%let path = "O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP";

/* 	result will be saved in folder	 */
libname HPV "O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP";


/** Call macros in follwoing sas program will do analysis for each data in the validation datasets **/
%include "O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\SAS code\Validation_HPV_macro.sas";


/** this is the Sars spike data spread sheets; */
%let file= HSL_VAL_HPV18L1L2VLP_Summary.xlsx;    


/****************************************************************************************************************/
/*****************************Do Analysis ***********************************************************************/
/****************************************************************************************************************/

/**************************************************LLOQ*********************************************************/


%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file", dat=LLOQ, sheet=LLOQ);

%lloq(dat=LLOQ);
options orientation=landscapes;
ODS RTF file="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HPV_18_L1L2VLP_LLOQ.rtf";

proc print data=Lloq_final1 noobs split='*' ;
  	 var HPV_type sample_id geomean Cv mean_result pct_err  ;
     label geomean='Geometric Mean(AU/mL)*' Cv='Coefficient of Variaton*' mean_result='concentration(AU/mL)*' pct_err='percent error*';
     title "Lower Limit of Quantitation for each sample";
	 run;

proc print data=lloq_sum noobs split='*' ;
  	 var HPV_type mean_conc median SD min max UCL LLOQ ;
     label mean_conc='Geometric Mean(AU/mL)*' SD='Standard Deviation*' UCL='Upper 95% CI*' LLOQ='LLOQ(AU/mL)*';
     title "Lower Limit of Quantitation";
	 run;
ODS RTF close;

/* remove some sampls and re-do analysis */

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HSL_VAL_HPV18L1L2VLP_Summary_1.xlsx", dat=LLOQ, sheet=LLOQ);

%lloq(dat=LLOQ);




/*************************************************ULOQ***********************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file", dat=ULOQ, sheet=ULOQ);


%lloq(dat=ULOQ);
ODS RTF file="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HPV_18_L1L2VLP_ULOQ_final.rtf";

proc print data=uloq_final1 noobs split='*' ;
  	 var HPV_type sample_id geomean Cv mean_result pct_err  ;
     label geomean='Geometric Mean(AU/mL)*' Cv='Coefficient of Variaton*' mean_result='concentration(AU/mL)*' pct_err='percent error*';
     title "Lower Limit of Quantitation for each sample";
	 run;

proc print data=uloq_sum noobs split='*' ;
  	 var HPV_type mean_conc median SD min max  ULOQ ;
     label mean_conc='Geometric Mean*' SD='Standard Deviation*'  ULOQ='ULOQ*';
     title "Upper Limit of Quantitation";
	 run;
ODS RTF close;


/* remove some sampls and re-do analysis */

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HSL_VAL_HPV18L1L2VLP_Summary_1.xlsx", dat=ULOQ, sheet=ULOQ);

%lloq(dat=ULOQ);

/**********************************************  Cutpoint  ********************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file",dat=CutPoint, sheet=CutPoint);

%cut(dat=CUTPOINT);

ODS RTF file="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HPV_18_L1L2VLP_cutpoint_final.rtf";
proc print data=Cutpoint_sum noobs split='*' ;
  	 var HPV_type mean median min max SD cutpoint;
     title "Cutpoint for Covid_19 NeucleoCapsid";
	 run;
ODS RTF close;


/********************************************  Precision ***********************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file", dat=PRECISION, sheet=PRECISION);

/********* calculate within day, between day and between analyst CVs**************/

%pres(dat=PRECISION);

ODS rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HPV_18_L1L2VLP_precision.rtf";
    proc print data=icc noobs split='*';
	 var Assay cvwd cvbd cvbt cv icc;
     label cvwd='Within Day CV*' cvbd='Between Day CV*' cvbt='Between Tech CV*' cv='Overall CV*' icc='ICC*';
     title "CVs and ICC for Neucleocapsid_precision";
	 run;
ods rtf close;




/****************************************************  accuracy **********************************************************************/

%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file", dat=ACCURACY, sheet=ACCURACY);

proc sql;
create table geo_analsyt_ACCURACY as
select unique HPV_type, Analyst, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV, (100-exp(mean(lg_acon))) as pct_err
				 /* analyst, day ,  Sample_id,*/
from ACCURACY
/* where lg_acon^=.  */
group by HPV_type,   Analyst, Sample_id;
quit;


proc sql;
create table geo_ACCURACY as
select unique HPV_type, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ 100*std(exp(lg_acon)) /mean(exp(lg_acon))  as CV, (100-exp(mean(lg_acon))) as pct_err
				 /* analyst, day ,  Sample_id,*/
from ACCURACY
/* where lg_acon^=.  */
group by HPV_type,    Sample_id;
quit;

ODS rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\HPV_18_L1L2VLP_accuracy.rtf";
    proc print data=geo_analsyt_ACCURACY noobs split='*';
	 var HPV_type analyst num_repplicate geomean cv pct_err;
     label geomean='Geometric mean(AU/mL)*' cv='%CV*' pct_err ='Percent Error*' cv='Overall CV*' icc='ICC*';
     title "HPV_18_L1L2VLP_accuracy by analyst";
	 run;

  proc print data=geo_ACCURACY noobs split='*';
	 var HPV_type num_repplicate geomean cv pct_err;
     label geomean='Geometric mean(AU/mL)*' cv='%CV*' pct_err ='Percent Error*' cv='Overall CV*' icc='ICC*';
     title "HPV_18_L1L2VLP_accuracy";
	 run;
ods rtf close;




/****************************************************************************************************************************************/
/**************************************************** Lot to Lot *************************************************************************/
/****************************************************************************************************************************************/

/*****************Antigen***********************/
%importdat(datfile="O:\HSL\HSL_COVID-19\Chunming Zhu\HPV\HPV_18_L1L2VLP\&file", dat=Stability, sheet=Stability);

%lot(dat=Stability);   /*    p=0.9761  */

option orientation =landscape;

ODS rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\Chunming_result\Sars_spike\Sars_spike_LOT_LOT_ANTIGEN.rtf";

proc print data=summ noobs split='*';
	 var HPV_type Num_exp mean_pct_err max_pct_err num_pass pass_pct;
	 label Num_exp='Antigen total*' mean_pct_err='Average percent of Error*' max_pct_err="Largest percent of Error*" num_pass='Antigen passing*' pass_pct ='Antigen with err <=25% *';;
     	 run;
ODS rtf CLOSE;


/*****************Conjugate***********************/
proc import datafile= 'O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\COVID19 Spike Validation Summary Final 14JUL21.xlsx'
out = Lot_c
dbms = xlsx
replace;
sheet = "Lot_to_Lot_Conjugate";
run;


%lot(dat=Lot_c);   /* p=0.8387   */

ODS rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\SARS\Chunming_result\Sars_spike\Sars_spike_LOT_LOT_COAGULANT.rtf";
proc print data=summ noobs split='*';
	 var HPV_type Num_exp mean_pct_err max_pct_err num_pass pass_pct;
	 label Num_exp='Coagulant total*' mean_pct_err='Average percent of Error*' max_pct_err="Largest percent of Error*" num_pass='Coagulant passing*' pass_pct ='Coagulant with err <=25% *';;
     	 run;
ODS rtf CLOSE;

