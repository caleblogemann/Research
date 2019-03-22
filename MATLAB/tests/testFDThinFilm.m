% for f = q
% fourth order operator
%lambda = 2*pi;
%exactSolution = @(x, t) exp(-1.0*lambda^4*t).*sin(lambda*x);
%exactSolution = @(x, t) 0.2*exp(-100.0*(x - 0.5).^2) + 0.1;
%forcingFunction = @(x) 8.0*exp(-100.0*(1.0 - 2.0*x).^2).*...
%    (2.0 + exp(25.0*(1.0 - 2.0*x).^2)).^2.*...
%    (exp(25.0*(1.0-2.0*x).^2).*(2203 - 18800*x + 58800*x.^2 - 80000*x.^3 + ...
%    40000*x.^4) + 2*(9253 - 77000*x + 237000*x.^2 - 320000*x.^3 + 160000*x.^4));

% time dependent exact solution
% exactSolution = @(x, t) 0.2*exp(-10*t).*exp(-100.0*(x - 0.5).^2) + 0.1;
% forcingFunction = @(x, t) exp(-40.0*(t+10.0*(x-0.5).^2)).*...
%     (288.0*exp(10.0*(t+10.0*(x-0.5).^2)).*...
%     (2301 - 19200*x + 59200*x.^2 - 80000*x.^3 + 40000*x.^4)+...
%     48.0*exp(20.0*(t+10.0*(x-0.5).^2)).*...
%     (4553 - 38200*x + 118200*x.^2 - 160000*x.^3 + 80000*x.^4)+...
%     2.0*exp(30.0*(t+10.0*(x-0.5).^2)).*...
%     (8811 - 75200*x + 235200*x.^2 - 320000*x.^3 + 160000*x.^4)+...
%     64*(9253 - 77000*x + 237000*x.^2 - 320000*x.^3 + 160000*x.^4));

exactSolution = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;
forcingFunction = @(x, t) 2.0*exp(-40.0*(t+30.0*(x-0.5).^2)).*...
    (648.0*exp(20.0*(t+30.0*(x-0.5).^2)).*...
    (14551 - 118200*x + 358200*x.^2 - 480000*x.^3 + 240000*x.^4)+...
    1296*exp(10.0*(t+30.0*(x-0.5).^2)).*...
    (21901 - 177600*x + 537600*x.^2 - 720000*x.^3 + 360000*x.^4)+...
    864.0*(29251 - 237000*x + 717000*x.^2 - 960000*x.^3 + 480000*x.^4)+...
    exp(30.0*(t+30.0*(x-0.5).^2)).*...
    (777707 - 6350400*x + 19310400*x.^2 - 25920000*x.^3 + 12960000*x.^4));


a = 0.0;
b = 1.0;

diffusivity = 1.0;
bc = 'periodic';

tFinal = 0.1;
error = [];
multiplier = 1;
timeStepArray = multiplier*[20, 40, 80, 160, 320, 640, 1280];
nCellsArray = [20, 40, 80, 160, 320, 640, 1280];
for i = 1:length(timeStepArray)
    nTimeSteps = timeStepArray(i);
    deltaT = tFinal/nTimeSteps;
    
    nCells = nCellsArray(i);
    deltaX = (b-a)/nCells;
    x = (a + deltaX/2):deltaX:(b-deltaX/2);
    forcingFunctionVector = @(t) forcingFunction(x',t);
    %forcingFunctionVector = @(t) zeros(nCells, 1);

    nBasisCpts = 2;
    q0 = projectQ(@(x) exactSolution(x, 0), nBasisCpts, nCells, a, b);
    qOnes = projectQ(@(x) ones(size(x)), nBasisCpts, nCells, a, b);
    qExact = @(t) projectQ(@(x) exactSolution(x, t+deltaT/2), 1, nCells, a, b);
    getAMatrix = @(q, t) getFDThinFilmMatrix(qExact(t), nCells, deltaX, diffusivity, bc);
    qFinal = IRK2(getAMatrix, q0, forcingFunctionVector, deltaT, tFinal);

    error = [error, norm(qFinal - exactSolution(x, tFinal)')/norm(exactSolution(x, tFinal)')];
    plot(x, qFinal, x, exactSolution(x, tFinal));
    disp(i);
    %pause()
end
deltaTArray = tFinal./timeStepArray;
deltaXArray = (b-a)./nCellsArray;
log(error(1:end-1)./error(2:end))./log(deltaTArray(1:end-1)./deltaTArray(2:end))
%log(error(1:end-1)./error(2:end))./log(deltaXArray(1:end-1)./deltaXArray(2:end))