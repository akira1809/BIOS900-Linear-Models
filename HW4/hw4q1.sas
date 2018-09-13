dm 'log; clear; output; clear;';

ods html close;
ods listing;

/*question1 */
proc iml;
	X1 = {1 0 -50 2500 ,
         1 0 0 0 ,
		1 0 50 2500,
		0 1 -50 2500,
		0 1 0 0,
		0 1 50 2500};
print X1;
rankX1=round(trace(ginv(X1)*X1));
  print X1 rankX1;
    X2 = {1 1 0 -1 1,
		  1 1 0 0 -2,
		  1 1 0 1 1,
		  1 0 1 -1 1,
		  1 0 1 0 -2,
          1 0 1 1 1};
	rankX2=round(trace(ginv(X2)*X2));
  print X2 rankX2;
  x0 = {1 0 -30 900};
  print x0;
  var_x0 = x0*inv(t(X1)*X1)*t(x0);
  print var_x0;
  coef = 4.302653*sqrt(1 + var_x0);
  print coef;
quit;

