#Q1, check the correctness of SVD for A

A= matrix(
  c(10, -5, 2, -11, 6, -8),
  nrow=3,
  ncol=2,
  byrow = TRUE)
U = matrix(
    c(1/sqrt(3), 1/sqrt(2), 1/sqrt(6),
      1/sqrt(3), -1/sqrt(2), 1/sqrt(6),
      1/sqrt(3), 0, -2/sqrt(6)
      ),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
V = matrix(
   c(0.6,0.8,-0.8,0.6),
   nrow=2,
   ncol=2,
   byrow=TRUE
)
Lam = matrix(
  c(10*sqrt(3),0,
    0,5*sqrt(2),
    0,0),
  nrow=3,
  ncol=2,
  byrow=TRUE
)
U%*%Lam%*%t(V)
