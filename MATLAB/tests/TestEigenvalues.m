a = 0;
b = 10;

qleft = 0.3323;
qright = 0.1;
discontinuity = 5;
qRiemann = @(x) qleft*double(x <= discontinuity) + qright*double(x > discontinuity);

amplitude = 0.5;
offset = 0.1;
mu = 5.0;
sigma = 1;
qGaussian = @(x) amplitude*1/(sqrt(2*pi*sigma^2))*exp(-(x - mu)^2/(2*sigma^2)) + offset;

amplitude = 0.1;
offset = 0.2;
wavenumber = 2;
qSin = @(x) amplitude*sin(2*pi*wavenumber*x/(b -a)) + offset;

deltaT = 0.4;
diffusivity = 1.0;
bc = 'extrapolation';
qFun = qRiemann;
nCells = 7;
deltaX = (b - a)/nCells;
nBasisCpts = 1;
q = projectQ(qFun, nBasisCpts, nCells, a, b);
qConstant = projectQ(qRiemann, nBasisCpts, nCells, a, b);
matrixFunctionThinFilm = @(num_elems) LDGThinFilmBE(num_elems, nBasisCpts, q, deltaT, a, b, diffusivity, bc);
matrixFunctionThinFilmPreconditioner = @(num_elems) LDGThinFilmBE(num_elems, nBasisCpts, qConstant, deltaT, a, b, diffusivity, bc);

M = matrixFunctionThinFilmPreconditioner(nCells);
A = matrixFunctionThinFilm(nCells);
[L, U] = lu(A);
[L, U] = ilu(sparse(M), struct('type', 'ilutp', 'droptol', 1e-5));
rhs = q;

%e = eig(full(A));
%plot(real(e), imag(e), 'o')
cA = cond(A);
cMIA = cond(M\A);
tol = 1e-8;
maxit = nCells*nBasisCpts;
[x, flag, relres, iter, resvec] = gmres(A, rhs, [], tol, maxit);
disp(iter(2))
[x, flag, relres, iter, resvec] = gmres(A, rhs, [], tol, maxit, M);
disp(iter(2))
[x, flag, relres, iter, resvec] = gmres(A, rhs, [], tol, maxit, L, U);
disp(iter(2))