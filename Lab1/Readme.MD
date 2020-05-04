**Ex 1.1.** Create a 10-dim row vector of all 1’s and then put to zero
its last 3 entries.

```
[1 1 1 1 1 1 1 1 1 1]→[1 1 1 1 1 1 1 0 0 0]
```
**Ex 1.2.** Create a 4×4 (uniformly) _random_ matrix with entries in
[0 _,_ 1] and then flip the elements on its _diagonal_.

**Ex 2.1.** Create a functionrvY = zeroTail(rvX)which has as
input an _n_ -dim ( _n_ ≥3) vectorrvXand as output a vectorrvY
equal torvXexcept for its last 3 entries which are set to 0.

**Ex 2.2.** Create a functionmY = flipDiag(mX)which has as input
an _n_ × _n_ matrixmXand returns a matrixmYequal tomXbut with a
flipped diagonal.

**Ex 2.3.** Create a functionbT = testSchur(rvP,rvQ)which has
as inputs two arbitrary polynomialsrvPandrvQ. This function:

**i)** plots the product of the two polynomials in the interval [− 10 _,_ 10],
**ii)** returns booleantrueif the latter product is Schur stable and
booleanfalseotherwise.

**Ex 3.1.** 

**i)** Create two (random) 500×2 measurement matricesmS1:= S 1 ,
mS2:= S 2 and stack them inmS:= S.
**ii)** Create a (random) 2×2 covariance matrixmP:= P , and two
(random) 500×500 covariance matricesmR1:= R 1 ,mR2:= R 2.
Use the latter matrices to build the matrixmR:= R.
**iii)** Generate a realization of the measurement vectorscvY1:= y 1 ,
cvY2:= y 2 and stack them incvY:= y.

**Ex 3.2.** Create the function
[cvE,mV] = centralMMSE(cvY,mP,mR,mS)
which has as input the generated realization cvY, the a pri-
ori covariancemP, the noise covariancemRand the measurement
matrixmS, and returns the best linear MMSE estimatecvEand the
variance of the estimation errormVusing the “standard” formula.

**Ex 3.3.** Create the function
[cvE,mV] = distribMMSE(cvY1,cvY2,mP,mR1,mR2,mS1,mS2)
which has as inputs the single-system measurement vectors,
noise covariance matrices and measurement matrices, and returns
the best linear MMSE estimate cvE and the variance of the
estimation errormVusing the “distributed” formula.
