dm'output;clear;log;clear;';


options  pageno=1  nodate  dsoptions=nomissopt;

title "senic data analysis";

ods html close;
ods listing;

libname data  "C:\akira\data";

/*import the data*/
proc import datafile = "C:\akira\data\senic.csv"
	dbms = dlm
	/*out = data.senic*/
	out = senic 
	replace;
	delimiter = ',';
	getnames = yes;
run;


proc print data = senic;
run;


/*to make scatter plots*/
ods graphics on;
ods html;
title 'senic data analysis';
proc corr data=senic nomiss plots=matrix(histogram);
   var length_of_stay age infection_risk available_facilities_and_service routine_chest_x_ray_ratio;
run;
ods graphics off;
ods html close;

proc iml;
/*use the senic data we imported*/
use work.senic;
/*read in the average length of stay into Y column*/
read all var {Length_of_stay} into Y;
/*read the 4 required variables into X matrix*/
read all var{age infection_risk available_facilities_and_service 
			routine_chest_x_ray_ratio} into X;

/*find the number of observations*/
d1 = nrow(X[, 1]);
print d1;

/*create column of 1s and insert into X*/
intercept = j(d1,1,1);
print intercept;
X_b = X;
X_a = intercept||X;

/*check dimension of design matrix*/
d_b = dimension(X_b);
d_a = dimension(X_a);
print d_b d_a;

/*create labels for matrices*/
cY = {"Y"};
cX_a = {"1" "X1" "X2" "X3" "X4"};
cX_b = {"X1" "X2" "X3" "X4"};
mattrib Y colname=cY X_a colname=cX_a X_b colname = cX_b;

/*optional steps, print to check Y and X*/
print Y;
print X_a;
print X_b;

/*OLS for beta*/
/*variables with _a is for part A model*/
/*variables with _b is for part B model*/
est_beta_a =inv(t(X_a)*X_a)*t(X_a)*Y;
est_beta_b = inv(t(X_b)*X_b)*t(X_b)*Y;
print est_beta_a;
print est_beta_b;

/*eigen values and eigen vectors of (X'X)^{-1}*/
a = inv(t(X_a)*X_a);
call eigen( e, u, a);
print,,, a,,, 'Eigenvalues and eigenvectors', e u;

/*compute the rank of X'X*/
a = t(X_a)*X_a;
b = t(X_b)*X_b;
ranka=round(trace(ginv(a)*a));
rankb=round(trace(ginv(b)*b));
call eigen(e_a, u_a, a);
call eigen(e_b, u_b, b);
print a ranka;
print e_a u_a;
print b rankb;
print e_b u_b;

/*compute the estimatd means est_Y*/
est_Y_a = X_a*est_beta_a;
est_Y_b = X_b*est_beta_b;
/*compute the residuals*/
res_a = Y - est_Y_a;
res_b = Y - est_Y_b;

print est_Y_a res_a;
print est_Y_b res_b;

print res_a res_b;


/*scatter plot of residual against estimated mean*/
ods html;
title "residual against estimated mean";
run Scatter(est_Y_a, res_a)
	/*add reference line*/
	other = "refline 0/axis = y"
	;
run Scatter(est_Y_b, res_b)
	/*add refernece line*/
	other = "refline 0/axis = y"
	;
ods html close;

/*export resiual to a SAS data set called work.residual*/
/*we make qq normal plot of residual outside of proc iml later*/
create residual var{res_a};
append;
close residual;

create residual_b var{res_b};
append;
close residual_b;

/*scatter qq plot of the residual*/
/*without using proc univariate*/
/*I need to understand quantile function*/

quit;


/*scattered normal qq plot for residual*/
title "Normal Quantile-Quantile Plot for Y residuals";
ods graphics on;
ods html;
proc univariate data=work.residual;
	qqplot res/odstitle=title;
run;
proc univariate data=work.residual_b;
	qqplot res_b/odstitle=title;
run;
ods graphics off;
ods html close;


