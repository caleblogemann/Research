global vcycle_depth weight num_relaxations tolerance max_iterations

weight = 2.0/3.0;
num_relaxations = 10;
tolerance = 1e-10;
max_iterations = 1000;

nBasisCpts = 1;
nCells = 1024;
vcycle_depth = 8;
a = 0;
b = 10;
deltaX = (b - a)/nCells;

qleft = 0.3323;
qright = 0.1;
discontinuity = 5.0;
qRiemann = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);

amplitude = 0.5;
offset = 0.1;
mu = 5.0;
sigma = 1;
qGaussian = @(x) amplitude*1/(sqrt(2*pi*sigma^2))*exp(-(x - mu)^2/(2*sigma^2)) + offset;

amplitude = 0.1;
offset = 0.2;
wavenumber = 2;
qSin = @(x) amplitude*sin(2*pi*wavenumber*x/(b -a)) + offset;

qFun = qRiemann;
q = projectQ(qFun, nBasisCpts, nCells, a, b);
rhs = q;
deltaT = 0.4;

diffusivity = 1.0;
matrixFunctionDiffusionFD = @(num_elems) FDDiffusionBE(num_elems, deltaT, a, b, diffusivity);
matrixFunctionDiffusion = @(num_elems) LDGDiffusionBE(num_elems, nBasisCpts, deltaT, a, b);
matrixFunctionHyperDiffusion = @(num_elems) LDGHyperDiffusionBE(num_elems, nBasisCpts, deltaT, a, b);
matrixFunctionThinFilm = @(num_elems) LDGThinFilmBE(num_elems, nBasisCpts, q, deltaT, a, b);
matrixFunction = matrixFunctionDiffusionFD;

[soln, num_iterations] = multigrid(q, rhs, matrixFunction);
exact_soln = getMatrix(matrixFunction(nCells)\getVector(rhs), nCells, nBasisCpts);
x = linspace(a, b, nCells);
plot(x, soln, x, exact_soln);
disp(num_iterations);
