dm'output;clear;log;clear;';


options  pageno=1  nodate  dsoptions=nomissopt;

title "senic data analysis";

ods html close;
ods listing;

/*Question 1*/
proc iml;
	a = {10 -5,
		  2 -11,
		  6 -8};
	call svd(u, q, v, a);
	print u, q, v;
quit;

dm'output;clear;log;clear;';

title "eruption data analysis";

ods html close;
ods listing;

libname data  "C:\akira\data";

/*Question 5*/
/*import the data*/
proc import datafile = "C:\akira\data\eruption.csv"
	dbms = dlm
/*	out = data.eruption*/
	out = work.eruption
	replace;
	delimiter = ',';
	getnames = yes;
run;

proc print data=eruption;
run;

proc iml;
use work.eruption;
read all var {y} into Y;
read all var{x} into X;

/*find the number of observations*/
d1 = nrow(X[, 1]);
print d1;

/*create column of 1s and insert into X*/
intercept = j(d1,1,1);
print intercept;
X_design = intercept||X;
print X_design;

/*check the dimension of design matrix*/
d = dimension(X_design);
print d;

/*create labels for matrices*/
cY = {"Y"};
cX = {"Intercept" "X"};
mattrib Y colname=cY X_design colname=cX;

/*check if the label works, this is optional*/
print Y;
print X_design;

/*part (a) OLS for beta*/
hat_beta = inv(t(X_design)*X_design)*t(X_design)*Y;
print hat_beta;
beta_0 = hat_beta[1, 1];
beta_1 = hat_beta[2, 1];
print beta_0;
print beta_1;

/*part (b) inference on beta_1*/

/*SSE, MSE, SSX and t*/
SSE = sum((y - X_design*hat_beta)##2);
print SSE;
MSE = SSE/(d1-2);
print MSE;
SSX = sum((x - mean(x))##2);
print SSX;

/*t statistic*/
t = beta_1/sqrt(MSE/SSX);
/*p value for alpha = 0.05*/
p = 2*(1 - probt(abs(t), d1-2));
print t;
print p;

/*confidence interval for beta_1*/
CI_lower = beta_1 - tinv(0.975, d1-2)*sqrt(MSE/SSX);
CI_upper = beta_1 + tinv(0.975, d1-2)*sqrt(MSE/SSX);
print CI_lower;
print CI_upper;

/*R-square: coefficient of determination*/
SSR = sum((X_design*hat_beta - mean(y))##2);
SST = sum((y - mean(y))##2);
R_square = SSR/SST;
print R_square;

quit;























