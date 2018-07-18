num_basis_cpts = 2;
num_elems = 10;
a = 0;
b = 1;
qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
deltaX = (b - a)/num_elems;
riemannFunction = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);
qRiemann = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
[LDG] = getLDGThinFilmMatrix(qRiemann, num_elems, num_basis_cpts, deltaX);
spy(LDG);
qRiemannVector = getQVector(qRiemann);
LDGqVector = LDG*qRiemannVector;

Aq = ThinFilmDiffusionOperator(qRiemann, deltaX, qRiemann, 'extrapolation');
AqVector = getQVector(Aq);
norm(AqVector - LDGqVector)

Dq = ThinFilmDiffusionDiagonalOperator(ones(num_elems, num_basis_cpts), deltaX, qRiemann, 'extrapolation');
DqVector = getQVector(Dq);
norm(diag(LDG) - DqVector)
