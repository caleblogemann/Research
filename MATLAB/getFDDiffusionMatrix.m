function [A] = FDDiffusionBE(nCells, deltaT, a, b, diffusivity)
    deltaX = (b - a)/nCells;
    M = diffusivity/deltaX^2;
    A = M
    R = getRMatrix(nCells, nBasisCpts, deltaX);
    Q = getQMatrixDiffusion(nCells, nBasisCpts, deltaX);

end
