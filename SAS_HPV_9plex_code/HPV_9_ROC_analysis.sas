/*********************************************************************/
/******validation for HPV9_*******************************************/
/*********************************************************************/


/*presicio data**/
/* use the geometric mean for *presicion* data each sample as expected value in LLOQ calculation*/

proc import datafile= 'O:/HSL/HSL_COVID-19/Chunming Zhu/Data analysis/HPV/HSL_VAL_HPV9_Summary_v5.xlsx'
out = pres
dbms = xlsx
replace;
sheet = "PRECISION";
run;

data pres;
set pres (rename= ( var10 = acon ) );
lg_acon = log(acon);
assay=hpv_type;

run;


proc sql;
create table geo_summary_pres as
select unique HPV_Type, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean
				 /* analyst, day ,*/
from pres
group by HPV_Type, Sample_id;
quit;


proc sort data=pres;
by assay sample_id  ;
run;




/****************************************************************************/
/****************************************cutpoint ***************************/
/****************************************************************************/
proc import datafile= 'O:/HSL/HSL_COVID-19/Chunming Zhu/Data analysis/HPV/HSL_VAL_HPV9_Summary_v5.xlsx'
out = cut
dbms = xlsx
replace;
sheet = "CUTPOINT";
run;

data cut;
set cut (rename= ( var10 = acon ) );
lg_acon = log(acon);
run;

proc sql;
create table geo_summary_cut as
select unique HPV_Type, Sample_id, count(*) as num_repplicate, exp(mean(lg_acon)) as geomean , /* Sample_id, */ exp(mean(lg_acon) +2*std(lg_acon)) as UCL
				 /* analyst, day ,  Sample_id,*/
from cut
group by HPV_Type, Sample_id;
quit;


/****************************************************************************/
/****************************  Lot  -  lot        ***************************/
/****************************************************************************/
proc import datafile= 'O:/HSL/HSL_COVID-19/Chunming Zhu/Data analysis/HPV/HSL_VAL_HPV9_Summary_v5.xlsx'
out = lot_B
dbms = xlsx
replace;
sheet = "Lot-to-Lot Bead Coupling";
run;

data lot_B;
set lot_B (rename= ( var10 = acon ) );
lg_acon = log(acon);
run;

proc sql;
create table geo_summary_lot_b as
select unique HPV_Type, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean
				 /* analyst, day ,*/
from lot_B
group by HPV_Type, Sample_id;
quit;


data lot_b_sum;
set geo_summary_lot_b ;
where sample_id in ("SS035" "SS038" "SS056" "SS094" "SS095" "SS107" "SS115" "SS113" "SS119" "HSL0001");
run;
/******************************************************************************************************************************/
proc import datafile= 'O:/HSL/HSL_COVID-19/Chunming Zhu/Data analysis/HPV/HSL_VAL_HPV9_Summary_v5.xlsx'
out = lot_PE
dbms = xlsx
replace;
sheet = "Lot-to-Lot Bead Conjugate (PE)";
run;

data lot_PE;
set lot_PE (rename= ( var10 = acon ) );
lg_acon = log(acon);
run;

proc sql;
create table geo_summary_lot_pe as
select unique HPV_Type, Sample_id,  count(*) as num_repplicate, exp(mean(lg_acon)) as geomean
				 /* analyst, day ,*/
from lot_PE
group by HPV_Type, Sample_id;
quit;


data lot_pe_sum;
set geo_summary_lot_pe;
where sample_id in ("SS035" "SS038" "SS056" "SS094" "SS095" "SS107" "SS115" "SS113" "SS119" "HSL0001");
run;

proc sql;
create table lot_sum_all as
select a.*, b.geomean as geomean2
from Lot_b_sum a left join Lot_pe_sum b on
a.HPV_Type = b.HPV_Type and a.sample_id = b.Sample_id
;
quit;
 

data lot_sum_all2;
set lot_sum_all;
geomean_f =exp(mean(log(geomean),log(geomean2)));
drop geomean geomean2;

run;


data lot_sum_all3;
set lot_sum_all2;
geomean=geomean_f;
drop geomean_f;
run;


/********************************************************************************/
/**          combine all ********************************************************/
/********************************************************************************/

data all ;
set Geo_summary_cut Geo_summary_pres lot_sum_all3;
run;


proc sort data=all;
by sample_id hpv_type;
run;




/***********************************************************************************/
/**********************Import sample information ***********************************/
/***********************************************************************************/

PROC IMPORT OUT= WORK.sample_infor 
            DATAFILE= "O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\HPV\sample_IDs.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc sort data=sample_infor ;
by sample_id;
run;

data all_merge;
merge sample_infor all;
by sample_id;
if sample_id="HSL00" then Response="Negative";
if Response="Negative" then resp =0;
	else resp =1;
if geomean=. then delete;
run;

proc sort data=all_merge;
by HPV_type resp geomean;
run;


proc freq data= all_merge;
table HPV_type* resp;
run;



/*******************************************************************************/
/****************   run model check ROC curve   ********************************/
/*******************************************************************************/




proc logistic data= all_merge plots=(ROC (id=prob));
by HPV_type;
class  resp(ref='0');
model resp = geomean;
output out=predata predicted=prob xbeta=beta;
run;


data respond;
set predata;
if . < prob <= 0.5 then resp_expected =0;
else if prob > 0.5 then resp_expected =1;
run;

ods rtf file="O:\HSL\HSL_COVID-19\Chunming Zhu\Data analysis\HPV\false_response.rtf";
proc freq data=respond;

table  HPV_type* resp*resp_expected/norow nocol nopercent nocum;
run;

ods rtf close;




data False_pred;
set predata ;
where (resp=1 and .< prob < 0.5) or (resp=0 and prob > 0.5);
run;

