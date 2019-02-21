function [q] = backwardEuler(getAMatrix, q0, forcingFunction, deltaT, tFinal)
    nTimeSteps = tFinal/deltaT;
    [nCells, ~] = size(q0);
    t = 0;
    q = q0;
    I = eye(nCells);
    for nT = 1:nTimeSteps
        A = getAMatrix(q, t);
        q = (I - deltaT*A)\(q + deltaT*forcingFunction(t + deltaT));
        t = t + deltaT;
    end
end
