function [LDG] = getLDGHyperDiffusionMatrix(nCells, nBasisCpts, deltaX, diffusivity, bc)
    R = getRMatrix(nCells, nBasisCpts, deltaX, bc);
    S = getSMatrix(nCells, nBasisCpts, deltaX);
    U = getUMatrixHyper(nCells, nBasisCpts, deltaX);
    Q = getQMatrix(nCells, nBasisCpts, deltaX, diffusivity);

    LDG = Q*U*S*R;
end
