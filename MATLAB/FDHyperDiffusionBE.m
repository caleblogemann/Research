function [A] = FDHyperDiffusionBE(nCells, deltaT, a, b, diffusivity, bc)
    deltaX = (b - a)/nCells;
    A = eye(nCells) - deltaT*getFDHyperDiffusionMatrix(nCells, diffusivity, deltaX, bc);
end
