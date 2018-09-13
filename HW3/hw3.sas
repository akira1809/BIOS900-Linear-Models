dm'output;clear;log;clear;';


options  pageno=1  nodate  dsoptions=nomissopt;

title "senic data analysis";

ods html close;
ods listing;

libname data  "C:\akira\data";

proc iml;

/*use the senic data we imported*/
use data.senic;
/*read in the average length of stay into Y column*/
read all var {Length_of_stay} into Y;
/*read the 4 required variables into X matrix*/
read all var{age infection_risk available_facilities_and_service 
			routine_chest_x_ray_ratio} into X;

/*set the number of parametsrs*/;
k = 4;

/*find the number of observations*/
n = nrow(X[, 1]);
print n;

/*create column of 1s and insert into X*/
intercept = j(n,1,1);
print intercept;
X_b = X;
X_a = intercept||X;

/*check dimension of design matrix*/
d_b = dimension(X_b);
d_a = dimension(X_a);
print d_b d_a;

/*compute beta hat*/
beta_hat = inv(t(X_a)*X_a)*t(X_a)*Y;
print beta_hat;

/*make the C matrix*/
c_1 = {0 1 0 0 0};
q_1 = 1;
c_2 = {0 1 0 0 0,
		0 0 0 1 0 };
q_2 = 2;
c_3 = {0 0 1 0 0, 
		0 0 0 1 0, 
	   0 0 0 0 1};
q_3 = 3;
c_4 = {0 1 0 -2 -2};
q_4 = 1;
c_5 = {0 -3 1 -3 -3};
q_5 = 1;
t = 0.25;

/*compute SSE and SSH*/

/*create identity matrix*/
I = i(n);
print I;
/*compute H matrix*/
H = X_a*inv(t(X_a)*X_a)*t(X_a);
/*compute SSE*/
SSE = t(Y)*(I - H)*Y;
print SSE;
/*compute SSH*/
SSH_1 = t(c_1*beta_hat)*inv((c_1*inv(t(X_a)*X_a)*t(c_1)))*(c_1*beta_hat);
SSH_2 = t(c_2*beta_hat)*inv((c_2*inv(t(X_a)*X_a)*t(c_2)))*(c_2*beta_hat);
SSH_3 = t(c_3*beta_hat)*inv((c_3*inv(t(X_a)*X_a)*t(c_3)))*(c_3*beta_hat);
SSH_4 = t(c_4*beta_hat)*inv((c_4*inv(t(X_a)*X_a)*t(c_4)))*(c_4*beta_hat);
SSH_5 = t(c_5*beta_hat - t)*inv((c_5*inv(t(X_a)*X_a)*t(c_5)))*(c_5*beta_hat - t);

/*compute F statistic and p-value*/:
F_1 = SSH_1/q_1/(SSE/(n - k - 1));
p_1 = 1 - CDF('F', F_1, q_1, n-k-1);
print F_1 p_1;

F_2 = SSH_2/q_2/(SSE/(n - k - 1));
p_2 = 1 - CDF('F', F_2, q_2, n-k-1);
print F_2 p_2;

F_3 = SSH_3/q_3/(SSE/(n - k - 1));
p_3 = 1 - CDF('F', F_3, q_3, n-k-1);
print F_3 p_3;

F_4 = SSH_4/q_4/(SSE/(n - k - 1));
p_4 = 1 - CDF('F', F_4, q_4, n-k-1);
print F_4 p_4;

F_5 = SSH_5/q_5/(SSE/(n - k - 1));
p_5 = 1 - CDF('F', F_5, q_5, n-k-1);
print F_5 p_5;

/*Simultaneous Tests*/

/*Bonferroni and */
a_1 = {0 1 0 0 0 };
a_2 = {0 0 1 0 0 };
a_3 = {0 0 0 1 0 };
a_4 = {0 0 0 0 1 };
/*compute F statistic and p values*/
F_critical_bon = finv((1 - 0.05/4), 1, n - k - 1);
F_a1 = t(a_1*beta_hat)*inv((a_1*inv(t(X_a)*X_a)*t(a_1)))*(a_1*beta_hat)/(SSE/(n - k - 1));
p_a1 = 1 - CDF('F', F_a1, 1, n - k - 1);
F_a2 = t(a_2*beta_hat)*inv((a_2*inv(t(X_a)*X_a)*t(a_2)))*(a_2*beta_hat)/(SSE/(n - k - 1));
p_a2 = 1 - CDF('F', F_a2, 1, n - k - 1);
F_a3 = t(a_3*beta_hat)*inv((a_3*inv(t(X_a)*X_a)*t(a_3)))*(a_3*beta_hat)/(SSE/(n - k - 1));
p_a3 = 1 - CDF('F', F_a3, 1, n - k - 1);
F_a4 = t(a_4*beta_hat)*inv((a_4*inv(t(X_a)*X_a)*t(a_4)))*(a_4*beta_hat)/(SSE/(n - k - 1));
p_a4 = 1 - CDF('F', F_a4, 1, n - k - 1);
print F_critical_bon;
print F_a1 p_a1 F_a2 p_a2 F_a3 p_a3 F_A4 p_a4;

/*Scheffe*/
/*compute (k + 1)F(alpha, k + 1,n - k - 1)*/
F_critical_s = (k + 1)*finv(1 - 0.05, k + 1, n- k - 1);
print F_a1 F_a2 F_a3 F_a4 F_critical_s;

/*individual confidence intervals for betas*/
g = inv(t(X_a)*X_a);
print g;
g11 = g[2, 2];
g22 = g[3, 3];
g33 = g[4, 4];
g44 = g[5, 5];
print g11 g22 g33 g44;
s = sqrt(SSE/(n - k - 1));
lower_1 = beta_hat[2] - tinv(0.975, n - k - 1)*s*sqrt(g11);
upper_1 = beta_hat[2] + tinv(0.975, n - k - 1)*s*sqrt(g11);
lower_2 = beta_hat[3] - tinv(0.975, n - k - 1)*s*sqrt(g22);
upper_2 = beta_hat[3] + tinv(0.975, n - k - 1)*s*sqrt(g22);
lower_3 = beta_hat[4] - tinv(0.975, n - k - 1)*s*sqrt(g33);
upper_3 = beta_hat[4] + tinv(0.975, n - k - 1)*s*sqrt(g33);
lower_4 = beta_hat[5] - tinv(0.975, n - k - 1)*s*sqrt(g44);
upper_4 = beta_hat[5] + tinv(0.975, n - k - 1)*s*sqrt(g44);
print lower_1 upper_1 lower_2 upper_2 lower_3 upper_3 lower_4 upper_4;

/*confidence interval for beta_1 - 2beta_3- 2beta_4*/
a = t({0 1 0 -2 -2});

lower_a = t(a)*beta_hat - tinv(0.975, n -  k - 1)*s*sqrt(t(a)*g*a);
upper_a = t(a)*beta_hat + tinv(0.975, n -  k - 1)*s*sqrt(t(a)*g*a);
print lower_a upper_a;

/*confidence interval for sigma^2*/
s_square = SSE/(n- k - 1);
lower_sigma = (n - k - 1)*s_square/cinv(0.975, n - k - 1);
upper_sigma = (n - k - 1)*s_square/cinv(0.025, n - k - 1);
print lower_sigma upper_sigma;

/*for part [C]*/

/*use only age X_1 as the predictor variable*/
read all var{age} into X_1;
X_1 = intercept||X_1;
beta_hat_1 = inv(t(X_1)*X_1)*t(X_1)*Y;
print beta_hat_1;
k_1 = 1;
H_1 = X_1*inv(t(X_1)*X_1)*t(X_1);
SSE_1 = t(Y)*(I - H_1)*Y;
print k_1;
s_1 = sqrt(SSE_1/(n - k_1 - 1));
/*given patient 67 years old*/
x_01 = {1 67};
/*compute prediction interval*/
lower_x_1 = x_01 * beta_hat_1 - tinv(0.975, n - k_1 - 1)*s_1*sqrt(1 + x_01*inv(t(X_1)*X_1)*t(x_01));
upper_x_1 = x_01 * beta_hat_1 + tinv(0.975, n - k_1 - 1)*s_1*sqrt(1 + x_01*inv(t(X_1)*X_1)*t(x_01));
print lower_x_1 upper_x_1;

/*for part [D]*/

/*compute likelihood ratio*/
LR = (1/(1 + q_4*F_4/(n - k - 1)))**(n/2);

test = -2*log(LR);
print test;
/*use chi square approximation for -2 log LR to obtain p value*/:
p_LR = 1 - probchi(-2*log(LR), 1);
print p_LR;

quit;

/*Question 3*/

/*create data*/
data Table7_5;
	input y x_1 x_2 x_3 @@;
	datalines;
	18.38	15.50	17.25	0.24	20.00	22.29	18.51	0.20
	11.50	12.36	11.13	0.12	25.00	31.84	 5.54	0.12
	52.50	83.90	 5.44	0.04	82.50	72.25	20.37	0.05
	25.00	27.14	31.20	0.27	30.67	40.41	 4.29	0.10
	12.00	12.42	 8.69	0.41	61.25	69.42	 6.63	0.04
	60.00	48.46	27.40	0.12	57.50	69.00	31.23	0.08
	31.00	26.09	28.50	0.21	60.00	62.83	29.98	0.17
	72.50	77.06	13.59	0.05	60.33	58.83	45.46	0.16
	49.75	59.48	35.90	0.32	 8.50	 9.00	 8.89	0.08
	36.50	20.64	23.81	0.24	60.00	81.40	 4.54	0.05
	16.25	18.92	29.62	0.72	50.00	50.32	21.36	0.19
	11.50	21.33	 1.53	0.10	35.00	46.85	 5.42	0.08
	75.00	65.94	22.10	0.09	31.56	38.68	14.55	0.17
	48.50	51.19	 7.59	0.13	77.50	59.42	49.86	0.13
	21.67	24.64	11.46	0.21	19.75	26.94	 2.48	0.10
	56.00	46.20	31.62	0.26	25.00	26.86	53.73	0.43
	40.00	20.00	40.18	0.56	56.67	62.52	15.89	0.05 
	;
run;

proc print data= Table7_5;
run;

proc iml;

/*use the senic data we imported*/
use work.Table7_5;
read all var {y} into Y;
read all var{x_1 x_2 x_3} into X;

/*set the number of parametsrs*/;
k = 3;

/*find the number of observations*/
n = nrow(X[, 1]);
print n;

/*create column of 1s and insert into X*/
intercept = j(n,1,1);
print intercept;
X = intercept||X;
print X;

/*check dimension of design matrix*/
d = dimension(X);
print d;

/*compute beta hat*/
beta_hat = inv(t(X)*X)*t(X)*Y;
print beta_hat;

/*create identity matrix*/
I = i(n);
print I;
/*compute H matrix*/
H = X*inv(t(X)*X)*t(X);
/*compute SSE*/
SSE = t(Y)*(I - H)*Y;
print SSE;
/*compute s*/
s = sqrt(SSE/(n - k - 1));

/*compute fitted values*/
y_hat = H*y;

/*compute the residuals*/

/*regular residuals*/
residual = (I - H)*Y;

/*studentized residuals*/
residual_stu = (inv(I - diag(H)))##(0.5)*residual/s;

/*deleted residuals*/
residual_del = inv(I - diag(H))*residual;
print residual_del;

/*external studentized residuals*/

/*compute SSEs without a single observation*/
/*refer to equation 9.32 from textbook*/
SSE_ext = SSE*j(n, 1, 1)  - inv(I - diag(H))*(residual##2);
/*compute s*/
s_ext = (SSE_ext/(n - k - 2))##(0.5);
/*compute esternal studentized residuals*/
residual_ext = (inv(I - diag(H))##(0.5))*(residual # (s_ext##(-1)));

/*scatter plot of residual against estimated mean*/
ods html;
title "studentized residuals against fitted value";
run Scatter(y_hat, residual_stu)
	/*add reference line*/
	other = "refline 0/axis = y"
	label={"fitted" "studentized res"}
	;
ods html close;
ods html;
title "deleted residuals against fitted value";
run Scatter(y_hat, residual_del)
	/*add refernece line*/
	other = "refline 0/axis = y"
	label={"fitted" "deleted res"}
	;
ods html close;
ods html;
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
leverage = diag(H)*j(n, 1, 1);
high_leverage = 2*(k + 1)/n;

/*compute cook distance*/
A = inv(I - diag(H))*diag(H);
D = A*(residual_stu##2)/(k + 1);

/*creating obervation number*/
obs= t(1:n);

/*creating table*/
table = obs||y||y_hat||residual||leverage||residual_stu||residual_ext||D;
/*create labels for our table*/
cTable = {"obs" "response" "fitted" "residual" "leverage(h_ii)" "r_i" 
		"t_i" "Cook"};
mattrib table colname=cTable; 
print table high_leverage PRESS SSE;

quit;

/*Question 4*/

proc iml;

/*create design matrix X*/

X = {1 1 0 0 , 1 1 0 0 ,
	 1 0 1 0 , 1 0 1 0 ,
	 1 0 1 0 , 1 0 1 0 ,
	 1 0 0 1,  1 0 0 1};

/*create matrix G*/
G = {0 0 0 0 , 0 0.5 0 0,
	 0 0 0.25 0 , 0 0 0 0.5};
print G;

/*Verify G is a generalized inverse of X'X*/
LHS = (t(X)*X)*G*(t(X)*X);
RHS = t(X)*X;
print LHS RHS;

/*get a solution for the MLE equation*/
A = G*t(X);
print A;

/*reparametrization*/
X_new ={1 0 0 , 1 0 0 ,
		0 1 0 , 0 1 0 ,
	    0 1 0 , 0 1 0 ,
		0 0 1, 0 0 1};
A_new = inv(t(X_new)*X_new)*t(X_new);
print A_new;
