function [A] = getFDHyperDiffusionMatrix(nCells, diffusivity, deltaX, bc)
    M = diffusivity/(deltaX^4);
    A = -6*M*eye(nCells) + 4*M*diag(ones(nCells-1,1),-1) + 4*M*diag(ones(nCells-1,1),1) - ...
         M*diag(ones(nCells-2, 1),2) - M*diag(ones(nCells-2,1),-2);
    if(strcmp(bc, 'extrapolation'))
        A(1,1) = -3*M;
        A(2,1) = 3*M;
        A(end-1,end) = 3*M;
        A(end, end)  = -3*M;
    elseif(strcmp(bc, 'periodic'))
        A(1,end-1) = -M;
        A(1,end)   = 4*M;
        A(2,end)   = -M;
        A(end-1,1) = -M;
        A(end,1)   = 4*M;
        A(end,2)   = -M;
    end
end
