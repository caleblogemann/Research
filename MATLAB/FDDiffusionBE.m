function [A] = FDDiffusionBE(nCells, deltaT, a, b, diffusivity)
    deltaX = (b - a)/nCells;
    M = diffusivity*deltaT/(deltaX^2);
    A = (2*M + 1)*eye(nCells) - M*diag(ones(nCells-1,1),-1) - M*diag(ones(nCells-1,1),1);
    A(1,1) = M + 1;
    A(end,end) = M + 1;

end
