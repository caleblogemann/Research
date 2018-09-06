function [A] = getFDDiffusionMatrix(nCells, diffusivity, deltaX, bc)
    M = diffusivity/(deltaX^2);
    A = (-2*M)*eye(nCells) + M*diag(ones(nCells-1,1),-1) + M*diag(ones(nCells-1,1),1);
    if(strcmp(bc, 'extrapolation'))
        A(1,1) = -M; 
        A(end,end) = -M;
    elseif(strcmp(bc, 'periodic'))
        A(end,1) = M;
        A(1,end) = M;
    end
end
