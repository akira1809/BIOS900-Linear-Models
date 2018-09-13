dm'output;clear;log;clear;';


options  pageno=1  nodate  dsoptions=nomissopt;

title "Question2: breaking strength data";

ods html close;
ods listing;
ods listing close;
ods html;

/*Question 2*/
proc iml;

/*	design matrix */
	W = {1 0 0 0 0 0 0 0 0 0 0 0 , 1 0 0 0 0 0 0 0 0 0 0 0 , 1 0 0 0 0 0 0 0 0 0 0 0,
		 0 1 0 0 0 0 0 0 0 0 0 0 , 0 1 0 0 0 0 0 0 0 0 0 0 , 0 1 0 0 0 0 0 0 0 0 0 0,
         0 0 1 0 0 0 0 0 0 0 0 0 , 0 0 1 0 0 0 0 0 0 0 0 0 , 
         0 0 0 1 0 0 0 0 0 0 0 0 , 0 0 0 1 0 0 0 0 0 0 0 0 , 0 0 0 1 0 0 0 0 0 0 0 0,
         0 0 0 0 1 0 0 0 0 0 0 0 , 0 0 0 0 1 0 0 0 0 0 0 0 , 0 0 0 0 1 0 0 0 0 0 0 0,
		 0 0 0 0 0 1 0 0 0 0 0 0 , 0 0 0 0 0 1 0 0 0 0 0 0 , 0 0 0 0 0 1 0 0 0 0 0 0, 0 0 0 0 0 1 0 0 0 0 0 0,
		 0 0 0 0 0 0 1 0 0 0 0 0 , 0 0 0 0 0 0 1 0 0 0 0 0 ,
         0 0 0 0 0 0 0 1 0 0 0 0 , 0 0 0 0 0 0 0 1 0 0 0 0 , 0 0 0 0 0 0 0 1 0 0 0 0,
         0 0 0 0 0 0 0 0 1 0 0 0 , 0 0 0 0 0 0 0 0 1 0 0 0 ,
         0 0 0 0 0 0 0 0 0 1 0 0 ,
         0 0 0 0 0 0 0 0 0 0 1 0 , 0 0 0 0 0 0 0 0 0 0 1 0 , 0 0 0 0 0 0 0 0 0 0 1 0,
         0 0 0 0 0 0 0 0 0 0 0 1 , 0 0 0 0 0 0 0 0 0 0 0 1 , 0 0 0 0 0 0 0 0 0 0 0 1 , 0 0 0 0 0 0 0 0 0 0 0 1};
	print W;
/*response y1*/
y = {21,27,19, 19, 19, 22, 19, 16, 23, 24, 23, 25, 23, 24, 23, 20, 24, 18, 19, 18, 28, 27, 25, 20, 24, 
	  28, 14, 16, 12, 23, 25, 22, 22};
/*create labels for matrices*/
cY = {"y"};
cW = {"mu_11" "mu_12" "mu_13" "mu_14" "mu_21" "mu_22" "mu_23" "mu_24" "mu_31" "mu_32" "mu_33" "mu_34"};
mattrib y colname=cY W colname=cW;
print W y;

/*contrasts matrix for testing main effect of cement(factor A), 
  main effect of aggregate(factor B), and interaction effect(A*B)*/
A = {2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1, 0 0 0 0 1 1 1 1 -1 -1 -1 -1};
print A;

B = {3 -1 -1 -1 3 -1 -1 -1 3 -1 -1 -1, 0 2 -1 -1 0 2 -1 -1 0 2 -1 -1,
	 0 0 1 -1 0 0 1 -1 0 0 1 -1};
print B;

C = {6 -2 -2 -2 -3 1 1 1 -3 1 1 1, 0 4 -2 -2 0 -2 1 1 0 -2 1 1,
     0 0 2 -2 0 0 -1 1 0 0 -1 1, 0 0 0 0 3 -1 -1 -1 -3 1 1 1,
     0 0 0 0 0 2 -1 -1 0 -2 1 1, 0 0 0 0 0 0 1 -1 0 0 -1 1};
print C;


/*compute OLS, H, SSE, s^2, SSA, SSB, SSAB*/
est_mu = inv(t(W)*W)*t(W)*y;
N = 33;
factor_A = 3;
factor_B = 4;
I = i(N);
H = W*inv(t(W)*W)*t(W);
SSE = t(y)*(I - H)*y;
v_E = N - factor_A * factor_B;
s_2 = SSE/(v_E);
print SSE s_2;

/*testing main effect of cement*/
SSA = t(A*estMu)*inv(A*inv(t(W)*W)*t(A))*A*estMu;
v_A = factor_A - 1;
/*F statistic and p value for testing main effect of cement*/
F_A = (SSA/v_A);
p_A = 1 - CDF('F', F_A, v_A, v_E);
print SSA F_A p_A;


/*testing main effect of aggreagate*/
SSB = t(B*estMu)*inv(B*inv(t(W)*W)*t(B))*B*estMu;
v_B = factor_B - 1;
/*F statistic and p value for testing main effect of cement*/
F_B = (SSB/v_B)/s_2;
p_B = 1 - CDF('F', F_B, v_B, v_E);
print SSB F_B p_B;

/*testing interaction between cement and aggregate*/
SSAB = t(C*estMu)*inv(C*inv(t(W)*W)*t(C))*C*estMu;
v_AB = (factor_A - 1)*(factor_B - 1);
/*F statistic and p value for testing main effect of cement*/
F_AB = (SSAB/v_AB)/s_2;
p_AB = 1 - CDF('F', F_AB, v_AB, v_E);
print SSAB F_AB p_AB;



/*for part (b)*/

/*compute K, A_new(A in the textbook), Z, Z1, K_star, mu_c, SSE_c, vE_c*/
j = {1 1 1 1 1 1 1 1 1 1 1 1 };
K = t(t(j)||t(A)||t(B));
A_new = t(t(K)||t(C));
print A_new;

rankA_new=round(trace(ginv(A_new)*A_new));
print rankA_new;

Z = W*inv(A_new);
Z1 = Z[, 1:6];
K_star = t(K)*inv(K*t(K));
est_mu_c = K_star*inv(t(Z1)*Z1)*t(Z1)*y;

SSE_c = t(y - W*est_mu_c)*(y - W*est_mu_c);
vE_c = v_E + 6;

/*test main effect of cement(A) and aggregate(B) for constrained model*/
/*output F statistic and p values*/

F_Ac = t(A*est_mu_c)*inv(A*K_star*inv(t(Z1)*Z1)*t(K_star)*t(A))*A*est_mu_c/v_A/SSE_c*vE_c;
p_Ac = 1 - CDF('F', F_Ac, v_A, vE_c);

F_Bc = t(B*est_mu_c)*inv(B*K_star*inv(t(Z1)*Z1)*t(K_star)*t(B))*B*est_mu_c/v_B/SSE_c*vE_c;
p_Bc = 1 - CDF('F', F_Bc, v_B, vE_c);

print SSE_c vE_c F_Ac p_Ac F_Bc p_Bc;


quit;

/*code for 2(c)*/
data q2;
	input cement aggregate strength @@;
	datalines;
	1 1 21		1 1 27		1 1 19		1 2 19		1 2 19		1 2 22
	1 3 19		1 3 16		1 4 23		1 4 24		1 4 23		2 1 25
	2 1 23		2 1 24		2 2 23		2 2 20		2 2 24		2 2 18
	2 3 19		2 3 18		2 4 28		2 4 27		2 4 25		3 1 20
	3 1 24		3 3 14		3 3 16		3 3 12		3 4 23		3 4 25
	3 4 22		3 4 22
	;
run;

proc print data=q2;
run;

proc glm data=q2;
	class cement aggregate;
	model strength = cement*aggregate/noint;
	contrast 'H_0 for cement' 
		cement*aggregate 2 0 2 2 -1 0 -1 -1 -1 -1 -1,
		cement*aggregate 0 0 0 0 1 0 1 1 -1 -1 -1;
	contrast 'H_0 for aggregate'
		cement*aggregate 3 -1 -1 -1 3 -1 -1 -1 0 0 0,
		cement*aggregate 0 2 -1 -1 0 2 -1 -1 0 0 0,
		cement*aggregate 0 0 1 -1 0 0 1 -1 0 0 0;
	contrast 'H_0 for cement-aggregate'
		cement*aggregate 1 -1 0 0 -1 1 0 0 0 0 0,
		cement*aggregate 0 1 -1 0 0 -1 1 0 0 0 0,
		cement*aggregate 0 1 0 -1 0 -1 0 1 0 0 0,
		cement*aggregate 1 0 -1 0 0 0 0 0 -1 1 0,
		cement*aggregate 1 0 0 -1 0 0 0 0 -1 0 1;
run;

proc glm data=q2;
	class cement aggregate;
	model strength = cement aggregate cement*aggregate/ss4;
run;

/*code for 2d*/
proc iml;
	/*	design matrix */
	W = {1 0 0 0 0 0 0 0 0 0 0 0 , 1 0 0 0 0 0 0 0 0 0 0 0 , 1 0 0 0 0 0 0 0 0 0 0 0 ,
		 0 1 0 0 0 0 0 0 0 0 0 0 , 0 1 0 0 0 0 0 0 0 0 0 0 , 0 1 0 0 0 0 0 0 0 0 0 0 ,
	     0 0 1 0 0 0 0 0 0 0 0 0 , 0 0 1 0 0 0 0 0 0 0 0 0 , 
         0 0 0 1 0 0 0 0 0 0 0 0 , 0 0 0 1 0 0 0 0 0 0 0 0 , 0 0 0 1 0 0 0 0 0 0 0 0 ,
         0 0 0 0 1 0 0 0 0 0 0 0 , 0 0 0 0 1 0 0 0 0 0 0 0 , 0 0 0 0 1 0 0 0 0 0 0 0 ,
         0 0 0 0 0 1 0 0 0 0 0 0 , 0 0 0 0 0 1 0 0 0 0 0 0 , 0 0 0 0 0 1 0 0 0 0 0 0 , 0 0 0 0 0 1 0 0 0 0 0 0 ,
		 0 0 0 0 0 0 1 0 0 0 0 0 , 0 0 0 0 0 0 1 0 0 0 0 0 ,
	     0 0 0 0 0 0 0 1 0 0 0 0 , 0 0 0 0 0 0 0 1 0 0 0 0 , 0 0 0 0 0 0 0 1 0 0 0 0 ,
         0 0 0 0 0 0 0 0 1 0 0 0 , 0 0 0 0 0 0 0 0 1 0 0 0 ,
		 0 0 0 0 0 0 0 0 0 1 0 0 , 0 0 0 0 0 0 0 0 0 1 0 0 , 0 0 0 0 0 0 0 0 0 1 0 0 ,
		 0 0 0 0 0 0 0 0 0 0 1 0 , 0 0 0 0 0 0 0 0 0 0 1 0 , 0 0 0 0 0 0 0 0 0 0 1 0 , 0 0 0 0 0 0 0 0 0 0 1 0,
		 0 0 0 0 0 0 0 0 0 0 0 0};
	print W;

	/*response y1*/
y = {21,27,19, 19, 19, 22, 19, 16, 23, 24, 23, 25, 23, 24, 23, 20, 24, 18, 19, 18, 28, 27, 25, 20, 24, 
	  14, 16, 12, 23, 25, 22, 22, 28};

/*create labels for matrices*/
cY = {"y"};
cW = {"mu_11" "mu_12" "mu_13" "mu_14" "mu_21" "mu_22" "mu_23" "mu_24" "mu_31" "mu_33" "mu_34" "mu_32" };
mattrib y colname=cY W colname=cW;
print W y;

/*contrasts matrix for testing main effect of cement(factor A), 
  main effect of aggregate(factor B), and interaction effect(A*B)*/
A = {2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1, 0 0 0 0 1 1 1 1 -1 -1 -1 -1};
print A;

B = {3 -1 -1 -1 3 -1 -1 -1 3 -1 -1 -1, 0 2 -1 -1 0 2 -1 -1 0 -1 -1 2,
      0 0 1 -1 0 0 1 -1 0 1 -1 0};
print B;

C = {6 -2 -2 -2 -3 1 1 1 -3 1 1 1 , 0 4 -2 -2 0 -2 1 1 0 1 1 -2, 
     0 0 2 -2 0 0 -1 1 0 -1 1 0, 0 0 0 0 3 -1 -1 -1 -3 1 1 1 ,
     0 0 0 0 0 2 -1 -1 0 1 1 -2, 0 0 0 0 0 0 1 -1 0 -1 1 0};
print C;

/*full reduced model approach*/
m = 1;
N = 33;
factor_A = 3;
factor_B = 4;
I = i(N);
H = W*ginv(t(W)*W)*t(W);
SSE_u = t(y)*(I - H)*y;
print SSE_u;

/*compute K, A_new(A in the textbook), Z, Z1, K_star, mu_c, SSE_c, vE_c*/
j = {1 1 1 1 1 1 1 1 1 1 1 1 };
K = t(t(j)||t(A)||t(B));
A_repar = t(t(K)||t(C));
print A_repar;

rankA_repar=round(trace(ginv(A_repar)*A_repar));
print rankA_repar;

K_star = t(K)*inv(K*t(K));
Z1 = W*K_star;
SSE_a = t(y)*(I - Z1*inv(t(Z1)*Z1)*t(Z1))*y;
F_int = (SSE_a - SSE_u)/(6 - 1)/SSE_u*(33-12);
p_int = 1 - CDF('F', F_int, 5, 21);
print  F_int p_int;

/*side condition approach*/
T = {1 -3 1 1 1 -3 1 1 -2 -2 -2 6};
/*equation 15.52*/
est_mu = inv(t(W)*W + t(T)*T)*t(W)*y;

SSE = t(y - W*est_mu)*(y - W*est_mu);
/*equation 15.53*/
cov_mu = inv(t(W)*W + t(T)*T)*t(W)*W*inv(t(W)*W + t(T)*T);
/*equation 15.54*/
F_int = t(C*est_mu)*ginv(C*cov_mu*t(C))*C*est_mu/(3*2 - 1)/SSE*(32 - 12 + 1);
p_int = 1 - CDF('F', F_int, 5, 21);
print SSE F_int p_int;

quit;

/*5(e)*/
proc iml;

	y = {21, 27, 19, 25, 23, 24, 20, 24, 19, 19, 22, 23, 20, 24, 18, 28,
		 19, 16, 19, 18, 14, 16, 12, 23, 24, 23, 28, 27, 25, 23, 25, 22, 22};
	N = nrow(y);

	/*cell means coding*/
	W1 ={1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    W2 = {0,0,0,0,0,0,0,0,1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	W3 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	W4 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	W = W1||W2||W3||W4;
	
    est_mu = inv(t(W)*W)*t(W)*y;
	SSB = t(est_mu)*t(W)*y - sum(y)##2/N;
	SSE = t(y)*y - t(est_mu)*t(W)*y;
	F = SSB/(4 - 1)/SSE*(N - 4);
	p = 1 - CDF('F', F, 3, N - 4);
	print SSE F p;

	/*the reference cell coding*/
	X1 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
		  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	X2 = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	X3 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    X4 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	      0, 1, 1, 1, 1, 1,1,1,1,1,1};
    X = X1||X2||X3||X4;
	
	est_mu = inv(t(X)*X)*t(X)*y;
	SSB = t(est_mu)*t(X)*y - sum(y)##2/N;
	SSE = t(y)*y - t(est_mu)*t(X)*y;
	F = SSB/(4 - 1)/SSE*(N - 4);
	p = 1 - CDF('F', F, 3, N - 4);
	print SSE F p;

	/*effect coding*/
	X1 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
		  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	X2 = {-1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	X3 = {-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 
		1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    X4 = {-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	      0, 1, 1, 1, 1, 1,1,1,1,1,1};
	X = X1||X2||X3||X4;
	
	est_mu = inv(t(X)*X)*t(X)*y;
	SSB = t(est_mu)*t(X)*y - sum(y)##2/N;
	SSE = t(y)*y - t(est_mu)*t(X)*y;
	F = SSB/(4 - 1)/SSE*(N - 4);
	p = 1 - CDF('F', F, 3, N - 4);
	print SSE F p;
	
quit;
