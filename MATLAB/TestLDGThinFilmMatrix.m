nBasisCpts = 2;
nCells = 10;
a = 0;
b = 1;
qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
deltaX = (b - a)/nCells;
riemannFunction = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);
qRiemann = projectQ(riemannFunction, nBasisCpts, nCells, a, b);
[LDG] = getLDGThinFilmMatrix(qRiemann, nCells, nBasisCpts, deltaX);
spy(LDG);
qRiemannVector = getQVector(qRiemann);
LDG*qRiemannVector

