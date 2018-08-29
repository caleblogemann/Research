function [LDG] = getLDGDiffusionMatrix(nCells, nBasisCpts, deltaX)
    R = getRMatrix(nCells, nBasisCpts, deltaX);
    Q = getQMatrixDiffusion(nCells, nBasisCpts, deltaX);

    LDG = Q*R;
end
