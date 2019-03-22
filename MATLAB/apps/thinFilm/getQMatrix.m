function [Q] = getQMatrix(nCells, nBasisCpts, deltaX, diffusivity)
    % Q_i = diffusivity/deltaX*(A - Phi(1, 1))U_i + diffusivity/deltaX*Phi(-1, 1)U_{i-1}
    % Q_i = H U_{i-1} + G U_i 
    % Q = [ H G 0 0 0 ]
    %     [ 0 H G 0 0 ]
    %     [ 0 0 H G 0 ]
    %     [ 0 0 0 H G ]
    %
    % add column to beginning of matrix for U_0
    % otherwise no boundary condition enforcement

    % one at top left
    sparseLowerDiagonal = speye(nCells, nCells+1);
    % one at bottom right
    sparseUpperDiagonal = sparse(1:nCells, 2:nCells+1, 1, nCells, nCells+1);
    
    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(1, 1, nBasisCpts);
    Phi2 = getPhiMatrix(-1, 1, nBasisCpts);
 
    G = diffusivity/deltaX*(A - Phi1);
    H = diffusivity/deltaX*Phi2;

    Q = kron(sparseLowerDiagonal, H) + kron(sparseUpperDiagonal, G);
end
