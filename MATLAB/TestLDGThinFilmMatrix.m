num_basis_cpts = 2;
num_elems = 50;
a = 0;
b = 1;
qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
diffusivity = 1.0;
bc = 'extrapolation';
deltaX = (b - a)/num_elems;
riemannFunction = @(x) qleft*double(x <= discontinuity) + qright*double(x > discontinuity);
height = 0.3;
mu = 0.5;
sigma = 0.1;
gaussianFunction = @(x) height*exp(-1/2*((x - mu)/sigma).^2);
qRiemann = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
[LDG] = getLDGThinFilmMatrix(qRiemann, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
spy(LDG);
qRiemannVector = getQVector(qRiemann);
LDGqVector = LDG*qRiemannVector;

Aq = ThinFilmDiffusionOperator(qRiemann, deltaX, qRiemann, diffusivity, 'extrapolation');
AqVector = getQVector(Aq);
norm(AqVector - LDGqVector)

Dq = ThinFilmDiffusionDiagonalOperator(ones(num_elems, num_basis_cpts), deltaX, qRiemann, diffusivity, 'extrapolation');
DqVector = getQVector(Dq);
norm(diag(LDG) - DqVector)
