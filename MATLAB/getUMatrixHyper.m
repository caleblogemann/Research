function [U] = getUMatrixHyper(nCells, nBasisCpts, deltaX)
    % U_i = -1/deltaX*(A + Phi(-1, -1))S_i + 1/deltaX*Phi(1, -1)S_{i+1}
    % U_i = E S_i + F S_{i+1}
    % U = [ E F 0 0 0 0 ]
    %     [ 0 E F 0 0 0 ]
    %     [ 0 0 E F 0 0 ]
    %     [ 0 0 0 E F 0 ]
    %     [ 0 0 0 0 E F ]
    %     Add first row to compute U_0 for Q
    %     Add last column to use S_{M+1}
    %     Added first column to use S_{0}
    %     U is (M+1) x (M+2) because added row and 2 columns

    % one at top left
    sparseLowerDiagonal = speye(nCells+1, nCells+2);
    % one at bottom right
    sparseUpperDiagonal = sparse(1:nCells+1, 2:nCells+2, 1, nCells+1, nCells+2);

    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = getPhiMatrix(1, -1, nBasisCpts);

    E = -1/deltaX*(A + Phi1);
    F = 1/deltaX*Phi2;

    U = kron(sparseLowerDiagonal, E) + kron(sparseUpperDiagonal, F);
end
