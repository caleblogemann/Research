% for f = q
% fourth order operator
%lambda = 2*pi;
%exactSolution = @(x, t) exp(-1.0*lambda^4*t).*sin(lambda*x);
exactSolution = @(x, t) 0.2*exp(-100.0*(x - 0.5).^2) + 0.1;
forcingFunction = @(x) 8.0*exp(-100.0*(1.0 - 2.0*x).^2).*...
    (2.0 + exp(25.0*(1.0 - 2.0*x).^2)).^2.*...
    (exp(25.0*(1.0-2.0*x)).^2.*(2203 - 18800*x + 58800*x.^2 - 80000*x.^3 + ...
    40000*x.^4) + 2*(9253 - 77000*x + 237000*x.^2 - 320000*x.^3 + 160000*x.^4));

a = 0.0;
b = 1.0;
nCells = 100;
deltaX = (b-a)/nCells;
x = (a + deltaX/2):deltaX:(b-deltaX/2);
forcingFunctionVector = @(t) forcingFunction(x');
%forcingFunctionVector = @(t) zeros(nCells, 1);
diffusivity = 1.0;
bc = 'periodic';

tFinal = 0.001;
error = [];
timeStepArray = [1, 2];
for nTimeSteps = timeStepArray
    
    deltaT = tFinal/nTimeSteps;

    nBasisCpts = 1;
    q0 = projectQ(@(x) exactSolution(x, 0), nBasisCpts, nCells, a, b);
    qOnes = projectQ(@(x) ones(size(x)), nBasisCpts, nCells, a, b);

    getAMatrix = @(q, t) getFDThinFilmMatrix(q, nCells, deltaX, diffusivity, bc);
    qFinal = backwardEuler(getAMatrix, q0, forcingFunctionVector, deltaT, tFinal);

    error = [error, norm(qFinal - exactSolution(x, tFinal)')];
    plot(x, qFinal, x, exactSolution(x, tFinal));
end
deltaTArray = tFinal./timeStepArray;
log(error(1:end-1)./error(2:end))./log(deltaTArray(1:end-1)./deltaTArray(2:end))