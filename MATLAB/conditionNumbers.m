qleft = 0.3323;
qright = 0.1;
discontinuity = 1;
qRiemann = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);
qSin = @(x) 0.1*sin(x) + 0.2;
mu = 1;
sigma = 0.5;
qGaussian = @(x) 1/(sigma*sqrt(2*pi))*exp(-1/2*((x - mu)/sigma).^2);
qFunctionArray = {qRiemann, qSin, qGaussian};

nCellsArray = [10, 20, 40, 80, 160];
nBasisCptsArray = [1, 2, 3];
conditionNumberMatrix = zeros(length(qFunctionArray), length(nCellsArray), length(nBasisCptsArray));
a = 0;
b = 2;
for qIndex = 1:length(qFunctionArray)
    qFunction = qFunctionArray{qIndex};
    for nCellsIndex = 1:length(nCellsArray)
        nCells = nCellsArray(nCellsIndex);
        deltaX = (b - a)/nCells;
        for nBasisCptsIndex = 1:length(nBasisCptsArray)
            nBasisCpts = nBasisCptsArray(nBasisCptsIndex);
            q = projectQ(qFunction, nBasisCpts, nCells, a, b);
            [LDG] = getLDGThinFilmMatrix(q, nCells, nBasisCpts, deltaX);
            conditionNumberMatrix(qIndex, nCellsIndex, nBasisCptsIndex) = cond(full(LDG));
        end
    end
end