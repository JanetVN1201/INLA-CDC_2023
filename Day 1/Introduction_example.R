library(INLA)
Q_10 = INLA:::inla.rw1(10)
diag(Q_10) <- diag(Q_10)+0.01
Q_10
Q_10_ns = matrix(Q_10, nrow = 10, byrow = T)
Q_10_ns
S_10 = SparseM::solve(Q_10)
S_10
S_10_ns = solve(Q_10_ns)
S_10_ns


##Create your own
i <- c(1,2,3,4,8)
j <- c(2:6)
x <- 1:5
A <- sparseMatrix(i, j, x = x)
A
summary(A)




