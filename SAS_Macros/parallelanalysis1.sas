
/* Either way, it seems simple. Here is what I did and it seems right so I am not clear why other macros I see are much more complicated. Please chime in if you see what I’m missing.

Randomly generate a set of random data with N variables and Y observations.
Keep the eigenvalues.
Repeat 500 times.
Combine the 500 datasets  (each will only have 1 record with N variables)
Find the 95th percentile

*/
%macro para(numvars,numreps) ;
%DO k = 1 %TO 500 ;
data A;
array nums {&numvars} a1- a&numvars ;
do i = 1 to &numreps;
do j = 1 to &numvars ;
nums{j} = rand(“Normal”) ;
if j < 2 then nums{j} = round(100*nums{j}) ;
else nums{j} = round(nums{j}) ;
end ;
drop i j ;
output;
end;

proc factor data= a outstat = a&k  noprint;
var a1 – a&numvars ;
data a&k ;
set  a&k  ;
if trim(_type_) = “EIGENVAL” ;

%END ;
%mend ;

%para(30,1000) ;

data all ;
set a1-a500 ;

proc univariate data= all noprint ;
var a1 – a30 ;
output out = eigvals pctlpts =  95 pctlpre = pa1 – pa30;

*** You don’t need the transpose but I just find it easier to read ;
proc transpose data= eigvals out=eigsig ;
Title “95th Percentile of Eigenvalues ” ;
proc print data = eigsig ;
run ;
/*
It runs fine and I have puzzled and puzzled over why a more complicated program would be necessary. I ran it 500 times with 1,000 observations and 30 variables and it took less than a minute on a remote desktop with 4GB RAM. Yes, I do see the possibility that if you had a much larger data set that you would want to optimize the speed in some way. Other than that, though, I can’t see why it needs to be any more complicated than this.

If you wanted to change the percentile, say, to 50, you would just change the 95 above. If you wanted to change the method from say, Principal Components Analysis (the default, with commonality of 1) to saying else, you could just do that in the PROC FACTOR step above.

The above assumes a normal distribution of your variables, but if that was not the case, you could change that in the RAND function above.

As I said, I am puzzled. Suggestions to my puzzlement welcome.
*/
