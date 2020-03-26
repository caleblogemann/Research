num_basis_cpts = 2;
num_cells = 10;
num_eqns = 1;
a = 0;
b = 2;
deltaX = (b - a)/num_cells;

diffusivity = 1.0;
bc = 'periodic';

qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
riemannFunction = @(x) qleft*double(x <= discontinuity) + ...
    qright*double(x > discontinuity);

height = 0.3;
mu = 0.5;
gaussian_sigma = 0.1;
gaussianFunction = @(x) height*exp(-1/2*((x - mu)/gaussian_sigma).^2);

sine_amplitude = 0.1;
sine_offset = 0.15;
sine_wavenumber = 1.0;
sine_wavespeed = 1.0;
sinFunction = @(x, t) sine_amplitude*sin(pi*sine_wavenumber*(x - ...
    sine_wavespeed*t)) + sine_offset;

quad_order = 5;
qSin = dog_math.L2Project(sinFunction, quad_order, num_cells, num_eqns, ...
    num_basis_cpts, a, b, 0.0);

q = qSin;
[LDG] = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, ...
    diffusivity, bc);

qVector = getQVector(q);
LDGqVector = LDG*qVector;

Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
AqVector = getQVector(Aq);
norm(AqVector - LDGqVector)

Dq = ThinFilmDiffusionDiagonalOperator(ones(num_elems, num_basis_cpts), ...
    deltaX, qRiemann, diffusivity, bc);
DqVector = getQVector(Dq);
norm(diag(LDG) - DqVector)
