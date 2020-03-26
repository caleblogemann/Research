function [A] = FDDiffusionBE(nCells, deltaT, a, b, diffusivity, bc)
    deltaX = (b - a)/nCells;
    A = eye(nCells) - deltaT*getFDDiffusionMatrix(nCells, diffusivity, deltaX, bc);
end
