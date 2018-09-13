dm 'log; clear; output; clear;';

ods html close;
ods listing;

data hw1q4;
input x1 x2 x3;
datalines;
3 10 4
2 11 5
6 15 8
6 9 9 
8 5 10
;
run;

proc means data = hw1q4;
run;

proc corr data = hw1q4 cov;
run;

proc iml;
/*standardize data: assume no missing values*/
start stdMat(x);
   mean = mean(x);                        /* means for columns */
   cx = x - mean;                     /* center x to mean zero */
   std = std(x);                 /* standard deviation estimate*/
   y = cx / std(x);                    /* scaling to std dev 1 */
   return( y );
finish stdMat;

/* Compute covariance: Assume no missing values*/
start covMat(x);
   n = nrow(x);                      /* number of observations */
   sum = x[+,];                         /* compute column sums */
   xpx = x`*x - sum`*sum/n;           /* compute sscp matrix   */
   /*notice that sscp is not covariance matrix, instead, sscp/(n- 1) is*/
   return(xpx/(n-1));
finish covMat;

/* Compute correlations: Assume no missing values  */
start corrMat(x);
   n = nrow(x);                      /* number of observations */
   sum = x[+,];                         /* compute column sums */
   xpx = x`*x - sum`*sum/n;           /* compute sscp matrix   */
   s = diag(1/sqrt(vecdiag(xpx)));           /* scaling matrix */
   corr = s*xpx*s;                       /* correlation matrix */
   return(corr);
finish corrMat;

/*the above code is not really necessary because we have cov() and corr() function
to deal with the same problem*/

use work.hw1q4;
read all var _NUM_ into x;
nm = {x1 x2 x3};
/*using self defined functions for cov, corr*/
/*std= stdMat(x);*/
/*cov = covMat(x);*/
/*corr = corrMat(x);*/

/*using SAS default functions for cov, corr*/
cov=cov(x);
corr = corr(x);
covinv = inv(cov(x));
id = cov*covinv;
stdx = mean(x)*inv(cov(x));
quadratic = mean(x)*inv(cov(x))*mean(x)`;
detcov = det(cov(x));
trace = trace(cov(x));
call eigen(eigenval, eigenvec, cov(x));
eigenval_square = eigval(cov(x)*cov(x));
sum_eig = sum(eigval(cov(x)));

/*compute the dot product of eigen vectors*/
eigvec_12 = sum(eigenvec[, 1]#eigenvec[, 2]);
eigvec_13 = sum(eigenvec[, 1]#eigenvec[, 3]);
eigvec_23 = sum(eigenvec[, 2]#eigenvec[, 3]);

/*solve homogeneous equation*/
a1 = {5 4 -2,
	  4 5 2,
	 -2 2 8	};
sol_1 = homogen(a1);
print sol_1; 

a2 = {-4 4 -2,
		4 -4 2,
		-2 2 -1};
sol_2 = homogen(a2);
print sol_2;

/* IML does not have a function that
   computes the rank of a matrix. You
   can use the following which makes 
   use of the generalized inverse 
   function  */

  a = {3 7 -3 6 4, 1 4 0 2 -5, 2 -1 1 4 -9, 0 -5 -3 0 -11, 0 -9 1 0 1};
  ranka=round(trace(ginv(a)*a));
  print a ranka;

  c = {7 -3 6 4, 4 0 2 -5 ,-1 1 4 -9 , -5 -3 0 -11};
  rankc=round(trace(ginv(c)*c));
  invc = inv(c);
  invctra = inv(c)`;
  print invc invctra;

  /*verify equation 2.58*/
  a = {3 7 -3 6 4, 1 4 0 2 -5, 2 -1 1 4 -9, 0 -5 -3 0 -11, 0 -9 1 0 1};
  aginv = {0 0 0 0 0, 
			-0.008333 0.2125 -0.09375 -0.022917 0,
			-0.108333 0.0125 0.15625 -0.172917 0,
			0.1 -0.175 0.1875 -0.0375 0,
			0.0333333 -0.1 0 -0.033333 0};
	checka = a*aginv*a;
	print checka;

/*verify 2.8c, question 7*/

/*part (i)*/

/*input A*/
a = {3 7 -3 6 4, 1 4 0 2 -5, 2 -1 1 4 -9, 0 -5 -3 0 -11, 0 -9 1 0 1};
/*getting generalized inverse*/
aginv = ginv(a);
/*getting rank of A*/
rankA = round(trace(ginv(a)*a));
/*getting rank of A inverse times A*/:
rankAinv_A = round(trace(ginv(ginv(a)*a)*(ginv(a)*a)));
/*getting rank of A times A inverse*/
rankA_Ainv = round(trace(ginv(a*ginv(a))*(a*ginv(a))));
/*print all three ranks to see if they are the same*/
print rankA rankAinv_A rankA_Ainv;

/*part (ii)*/
A_trans = a`;
/*computing transpose of inverse*/
invA_trans = ginv(a)`;
/*computing inverse of transpose*/
tranA_inv = ginv(a`);
/*print both to see if they are the same*/
print invA_trans;
print tranA_inv;

/*part (iii)*/

/*computing a transpose*/
Atran = a`;

/*computing A(A'A)^-A'A*/
a1 = a*ginv((a`*a))*a`*a;

/*computing (A'A)(A'A)^-A'*/
a2 = (a`*a)*ginv(a`*a)*a`;

/*print to see if they match*/
print a;
print a1;
print Atran;
print a2;

/*part (iv)*/

/*compute generalized inverse of A*/
Ainv = ginv(a);

/*compute (A'A)^-A'*/
Ainv_2 = ginv(Atran*a)*Atran;

print Ainv;
print Ainv_2;

/*part (v)*/

/*compute A(A'A)^-A' and its transpose*/
a3 = a*ginv(Atran*a)*Atran;
a4 = (a*ginv(Atran*a)*Atran)`;

/*compute the ran of A(A'A)^-A'*/
ranka3 = round(trace(ginv(a*ginv(Atran*a)*Atran)*a*ginv(Atran*a)*Atran));

print a3;
print a4;
print ranka3;

/*question 7 part C*/

/*part (i)*/

/*compute the trace of A'A and AA'*/
trace1 = trace(Atran*A);
trace2 = trace(A*Atran);

print trace1;
print trace2;

/*part (ii)*/

/*creating orthogonal matrix C*/
call comport(q,r,p,piv,lindep,A);
C = q;
/*check to see if C is orthogonal*/
/*annotate this in the formal code*/
identity = C*C`;
print C;
print identity;

/*compute the trace of C'AC and A*/

a5= C`*a*C;

trace3 = trace(a5);
trace4 = trace(a);

print trace3;
print trace4;

/*part (iii)*/

/*compute the trace for A^-A and AA^-*/
a6 = a*ginv(a);
a7 = ginv(a)*a;

/*compute their traces*/
trace6 = trace(a6);
trace7 = trace(a7);

print trace6;
print trace7;

/*part (iv)*/

/*compute sum of eigen values for A*/
trace8 = sum(eigval(a));

/*compute trace of A*/
traceA = trace(a);

print trace8;
print traceA;


print x;
print cov[rowname=nm colname=nm label= "Covariance Matrix"];
print corr[rowname=nm colname=nm label="Correlation Matrix"];
print covinv[rowname=nm colname=nm label="Covariance Inverse Matrix"];
print id[rowname=nm colname=nm label="Cov multiply invers cov"];
print std[colname = nm label = "Standardized Data"];
print stdx;
print quadratic;
print detcov;
print trace;
print 'Eigenvalues and eigenvectors of S', eigenval, eigenvec;
print 'eigen values of S^2', eigenval_square;
print 'sum of eigenvalues of S', sum_eig;
print 'pairwise dot product of eigenvectors:', eigvec_12, eigvec_13,eigvec_23; 
quit;

