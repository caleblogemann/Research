function [q] = IRK2(getAMatrix, q0, forcingFunction, deltaT, tFinal)
    nTimeSteps = tFinal/deltaT;
    [nCells, ~] = size(q0);
    t = 0;
    q = q0(:,1);
    I = eye(nCells);
    for nT = 1:nTimeSteps
        A = getAMatrix(q, t);
        qstar = (I - 0.25*deltaT*A)\((I + 0.25*deltaT*A)*q + 0.25*deltaT*(forcingFunction(t) + forcingFunction(t + 0.5*deltaT)));
        q = (3*I - deltaT*A)\(4*qstar - q + deltaT*forcingFunction(t + deltaT));
        t = t + deltaT;
    end
end
