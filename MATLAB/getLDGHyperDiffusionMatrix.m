function [LDG] = getLDGHyperDiffusionMatrix(nCells, nBasisCpts, deltaX)
    %qVector = getQVector(q);
    R = getRMatrix(nCells, nBasisCpts, deltaX);
    %rVector = R*qVector;
    %disp(rVector);
    S = getSMatrix(nCells, nBasisCpts, deltaX);
    %sVector = S*rVector;
    %disp(sVector);
    U = getUMatrixHyper(nCells, nBasisCpts, deltaX);
    %uVector = U*sVector;
    %disp(uVector);
    Q = getQMatrix(nCells, nBasisCpts, deltaX);

    LDG = Q*U*S*R;
end
