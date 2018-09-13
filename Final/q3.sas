dm'output;clear;log;clear;';


options  pageno=1  nodate  dsoptions=nomissopt;

title "senic data analysis";

ods html close;
ods listing;

/*Question 3*/
proc iml;
/*	design matrix */
	X = {1 162 23 3, 1 162 23 8, 1 162 30 5, 1 162 30 8,
		 1 172 25 5, 1 172 25 8, 1 172 30 5, 1 172 30 8,
         1 167 27.5 6.5, 1 177 27.5 6.5, 1 157 27.5 6.5, 1 167 32.5 6.5,
         1 167 22.5 6.5, 1 167 27.5 9.5, 1 167 27.5 3.5, 1 177 20 6.5,
		 1 177 20 6.5, 1 160 34 7.5, 1 160 34 7.5};
/*response y1*/
y1 = {41.5, 33.8, 27.7, 21.7, 19.9, 15.0, 12.2, 4.3, 19.3, 6.4, 37.6, 
			18.0, 26.3, 9.9, 25.0, 14.1, 15.2, 15.9, 19.6};
/*create labels for matrices*/
cY = {"y1"};
cX = {"1" "x1" "x2" "x3" };
mattrib y1 colname=cY X colname=cX;

print X y1;

rankX=round(trace(ginv(X)*X));

print rankX;

/*find dimension of design matrix X*/
dX =dimension(X);
print dX;

n = 19;
k = 3;
/*compute I,  H, SSE, hat_beta, and s^2*/
I = i(n);
H = X*inv(t(X)*X)*t(X);
SSE = t(y1)*(I - H)*y1;
s_2 = SSE/(n - k - 1);
hat_beta = inv(t(X)*X)*t(X)*y1;

print hat_beta s_2;

/*beta covariance matrix estimate*/
cov_beta_hat = s_2*inv(t(X)*X);

print cov_beta_hat;

/*compute SSR*/
SSR = t(hat_beta)*t(X)*y1 - sum(y1)##2/n;
print SSR;

/*R^2 and R^2_a*/
R_2 = SSR/(SSE + SSR);
R_2_adj = ((n - 1)*R_2 - k)/(n - k - 1);

print R_2 R_2_adj;

/*create second order terms*/

/*square terms*/
x1_2 = X[, 2]##2;
x2_2 = X[, 3]##2;
x3_2 = X[, 4]##2;
print x1_2 x2_2 x3_2;
/*mixed terms*/
x12 = X[, 2]#X[, 3];
x13 = X[, 2]#X[, 4];
x23=  X[, 3]#X[, 4];

print x1_2 x2_2 x3_2 x12 x13 x23;

X_new = X||x1_2||x2_2||x3_2||x12||x13||x23;

print X_new;

rankx_new = round(trace(ginv(X_new)*X_new));

print rankx_new;

k_new = 9;
H_new = X_new*inv(t(X_new)*X_new)*t(X_new);
SSE_new = t(y1)*(I - H_new)*y1;
s_2_new = SSE_new/(n - k_new - 1);
s_new = s_2_new##(0.5);
hat_beta_new = inv(t(X_new)*X_new)*t(X_new)*y1;
SSR_new = t(hat_beta_new)*t(X_new)*y1 - sum(y1)##2/n;
R_2_new = SSR_new/(SSE_new + SSR_new);
R_2_adj_new = ((n - 1)*R_2_new - k_new)/(n - k_new - 1);

print SSE_new s_2_new hat_beta_new SSR_new R_2_new R_2_adj_new;

/*residual, outlier and influential observation analysis*/

/*compute fitted values*/
y_hat = H_new*y1;

/*compute the residuals*/

/*regular residuals*/
residual = (I - H_new)*y1;

/*studentized residuals*/
residual_stu = (inv(I - diag(H_new)))##(0.5)*residual/s_new;

/*deleted residuals*/
residual_del = inv(I - diag(H_new))*residual;
print residual_del;

/*external studentized residuals*/

/*compute SSEs without a single observation*/
/*refer to equation 9.32 from textbook*/
SSE_ext = SSE_new*j(n, 1, 1)  - inv(I - diag(H_new))*(residual##2);
/*compute s*/
s_ext = (SSE_ext/(n - k_new - 2))##(0.5);
/*compute esternal studentized residuals*/
residual_ext = (inv(I - diag(H_new))##(0.5))*(residual # (s_ext##(-1)));

/*scatter plot of residual against estimated mean*/
ods html;
ods listing close;
title "studentized residuals against fitted value";
run Scatter(y_hat, residual_stu)
	/*add reference line*/
	other = "refline 0/axis = y"
	label={"fitted" "studentized res"}
	;
title "deleted residuals against fitted value";
run Scatter(y_hat, residual_del)
	/*add refernece line*/
	other = "refline 0/axis = y"
	label={"fitted" "deleted res"}
	;
title "ordinary residual against deleted residual";
run Scatter(residual, residual_del)
	/*add refernece line*/
	other = "refline 0/axis = y"
	label={"ordinary res" "deleted res"}
	;
ods html close;

ods listing;

/*comopute PRESS(prediction sum of square)*/
PRESS = t(residual_del)*residual_del;
print PRESS;

/*compute leverage*/
leverage = diag(H_new)*j(n, 1, 1);
high_leverage = 2*(k_new + 1)/n;

print high_leverage;

/*compute cook distance*/
A = inv(I - diag(H_new))*diag(H_new);
D = A*(residual_stu##2)/(k_new + 1);

/*creating obervation number*/
obs= t(1:n);

/*creating table*/
table = obs||y1||y_hat||residual||leverage||residual_stu||residual_ext||D;
/*create labels for our table*/
cTable = {"obs" "response" "fitted" "residual" "leverage(h_ii)" "r_i" 
		"t_i" "Cook"};
mattrib table colname=cTable; 
print table high_leverage PRESS SSE_new;

/*computing confidence intervals*/

/*individual confidence intervals for betas*/
g = inv(t(X_new)*X_new);
print g;
g11 = g[2, 2];
g22 = g[3, 3];
g33 = g[4, 4];
g44 = g[5, 5];
g55 = g[6, 6];
g66 = g[7, 7];
g77 = g[8, 8];
g88 = g[9, 9];
g99 = g[10, 10];
lower_1 = hat_beta_new[2] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g11);
upper_1 = hat_beta_new[2] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g11);
lower_2 = hat_beta_new[3] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g22);
upper_2 = hat_beta_new[3] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g22);
lower_3 = hat_beta_new[4] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g33);
upper_3 = hat_beta_new[4] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g33);
lower_4 = hat_beta_new[5] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g44);
upper_4 = hat_beta_new[5] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g44);
lower_5 = hat_beta_new[6] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g55);
upper_5 = hat_beta_new[6] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g55);
lower_6 = hat_beta_new[7] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g66);
upper_6 = hat_beta_new[7] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g66);
lower_7 = hat_beta_new[8] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g77);
upper_7 = hat_beta_new[8] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g77);
lower_8 = hat_beta_new[9] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g88);
upper_8 = hat_beta_new[9] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g88);
lower_9 = hat_beta_new[10] - tinv(0.975, n - k_new - 1)*s_new*sqrt(g99);
upper_9 = hat_beta_new[10] + tinv(0.975, n - k_new - 1)*s_new*sqrt(g99);
print lower_1 upper_1 lower_2 upper_2 lower_3 upper_3 lower_4 upper_4 lower_5 upper_5;
print lower_6 upper_6 lower_7 upper_7 lower_8 upper_8 lower_9 upper_9;

/*confidence interval for sigma^2*/
lower_sigma = (n - k_new - 1)*s_2_new/cinv(0.975, n - k_new - 1);
upper_sigma = (n - k_new - 1)*s_2_new/cinv(0.025, n - k_new - 1);
print lower_sigma upper_sigma;

quit;

/*Verify with regular SAS code*/
data ex7_54;
	input y1 x1 x2 x3 @@;
	datalines;
	41.5 162 23 3		33.8 162 23 8 		27.7 162 30 5		21.7 162 30 	8
	19.9 172 25 5		15.0 172 25 8		12.2 172 30 5		4.3  172 30 	8
	19.3 167 27.5 6.5	6.4  177 27.5 6.5	37.6 157 27.5 6.5	18.0 167 32.5 6.5
	26.3 167 22.5 6.5	9.9  167 27.5 9.5	25.0 167 27.5 3.5	14.1 177 20   6.5	
	15.2 177 20   6.5 	15.9 160   34 7.5   19.6 160  34  7.5
	;
run;

proc print data=ex7_54;

run;

proc reg data=ex7_54 outest = est covout;
	model y1 = x1-x3;
run;

proc print data=est;
run;

data ex7_54_new;
	set ex7_54;
	x4 = x1*x1;
	x5 = x2*x2;
	x6 = x3*x3;
	x7 = x1*x2;
	x8 = x1*x3;
	x9 = x2*x3;
run;

proc print data=ex7_54_new;
run;

proc reg data=ex7_54_new;
	model y1 = x1-x9/clb;
	output out = check
		   cookd = cook
		   PRESS = press;
run;

proc print data=check;
run;
