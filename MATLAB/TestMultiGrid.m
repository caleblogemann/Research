global vcycle_depth weight num_relaxations tolerance max_iterations

a = 0;
b = 10;

qleft = 0.3323;
qright = 0.1;
discontinuity = 5;
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

weight = 2.0/3.0;
num_relaxations = 10;
tolerance = 1e-6;
max_iterations = 1000;

deltaT = 0.4;
diffusivity = .01;
bc = 'extrapolation';


ndoublings = 4;
num_iterations_array = zeros(ndoublings,1);
initial_nCells = 64;
initial_vcycle_depth = 4;

nBasisCpts = 1;
qFun = qRiemann;

for i = 1:ndoublings
    nCells = initial_nCells*2^(i-1);
    deltaX = (b - a)/nCells;
    vcycle_depth = initial_vcycle_depth + (i-1);
    q = projectQ(qFun, nBasisCpts, nCells, a, b);
    rhs = q;
    matrixFunctionDiffusionFD = @(num_elems) FDDiffusionBE(num_elems, deltaT, a, b, diffusivity, bc);
    matrixFunctionDiffusion = @(num_elems) LDGDiffusionBE(num_elems, nBasisCpts, deltaT, a, b, diffusivity, bc);
    matrixFunctionHyperDiffusionFD = @(num_elems) FDHyperDiffusionBE(num_elems, deltaT, a, b, diffusivity, bc);
    matrixFunctionHyperDiffusion = @(num_elems) LDGHyperDiffusionBE(num_elems, nBasisCpts, deltaT, a, b, diffusivity, bc);
    matrixFunctionThinFilm = @(num_elems) LDGThinFilmBE(num_elems, nBasisCpts, q, deltaT, a, b, diffusivity, bc);
    matrixFunction = matrixFunctionHyperDiffusion;
    [soln, num_iterations_array(i)] = multigrid(q, rhs, matrixFunction);
    exact_soln = getMatrix(matrixFunction(nCells)\getVector(rhs), nCells, nBasisCpts);
    x = linspace(a, b, nCells);
    plot(x, soln, x, exact_soln);
    disp(num_iterations_array(i));
    disp(norm(soln - exact_soln));
end
