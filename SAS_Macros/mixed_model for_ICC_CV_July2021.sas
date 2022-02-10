
*    %let dat=pres;
/*
ODS output CovParms=cvtable;
PROC MIXED METHOD=REML data=&dat ;
*by assay;
CLASS date sid analyst;
MODEL AI50 =  /ddfm=satterth; 
RANDOM intercept /SUBJECT=Sid;
RANDOM date analyst /SUBJECT=Sid;
run;


data CV;
  set cvtable;
  if CovParm='Residual' then e=estimate;
  if CovParm='Intercept' then s=estimate;
  if CovParm='Date' then d=estimate;
  if CovParm='Analyst' then t=estimate;
  format s best32. e best32. d best32. t best32.;
run;
*/



%macro ICC_CV (dat);

ODS output CovParms=cvtable;
PROC MIXED METHOD=REML data=&dat ;
by assay;
CLASS day sample_id Analyst;
MODEL lg_acon =  /ddfm=satterth; 
RANDOM intercept  /SUBJECT=Sample_id;
RANDOM intercept /SUBJECT=day(analyst*sample_id) ;

repeated / subject=analyst(sample_id*day)  type=cs;
run;


data CV;
  set cvtable;
  if CovParm='Residual' then e=estimate;
  if Subject='Sample_ID' then s=estimate;
  if Subject='Day(Sample_*Analyst)' then d=estimate;
  if Subject='Analyst(Day*Sample_)' then t=estimate;
  format s best32. e best32. d best32. t best32.;
run;

proc univariate noprint data=CV; 
by assay;
var s d e t;
output out=CV2 sum=s d e t;
run;

data icc; 
    set cv2;
    cvwd=100*(abs(e)**.5);
    cvbd=100*(abs(d)**.5);
    cvbt=100*(abs(t)**.5);
    cv=100*(abs(e+d+t)**.5);
    icc=100*s/abs(s+d+e+t);
run;


*ODS HTML body="O:\HPV Assays Validation_Qualification Reports\HPV16 ELISA\SAS code\CVs and ICC for Overall reproducibility test.html";
     proc print data=icc noobs split='*';
	 var cvwd cvbd cvbt cv icc;
     label cvwd='Within Day CV*' cvbd='Between Day CV*' cvbt='Between Tech CV*' cv='Overall CV*' icc='ICC*';
     title "CVs and ICC for Overall reproducibility test";
	 run;
 *ODS HTML Close;
%mend;
