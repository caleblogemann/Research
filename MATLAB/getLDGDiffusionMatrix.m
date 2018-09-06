function [LDG] = getLDGDiffusionMatrix(nCells, nBasisCpts, deltaX, diffusivity, bc)
    R = getRMatrixDiffusion(nCells, nBasisCpts, deltaX, bc);
    Q = getQMatrixDiffusion(nCells, nBasisCpts, deltaX, diffusivity);

    LDG = Q*R;
end
