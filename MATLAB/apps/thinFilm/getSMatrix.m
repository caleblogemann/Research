function [S] = getSMatrix(nCells, nBasisCpts, deltaX)
    % S_i = - 1/deltaX*Phi(-1, 1)R_{i-1} - 1/deltaX*(A - Phi(1, 1))R_i 
    % S_i = C R_{i-1} + D R_i 
    % S = [ C D 0 0 0 0 0 ]
    %     [ 0 C D 0 0 0 0 ]
    %     [ 0 0 C D 0 0 0 ]
    %     [ 0 0 0 C D 0 0 ]
    %     [ 0 0 0 0 C D 0 ]
    %     [ 0 0 0 0 0 C D ]
    %
    % Add first row to compute S_0 and last row to compute S_{M+1} - needed for U
    % Added 2 first columns for R_{-1} and R_{0} needed for S_0 and S_1
    % Added last colum for R_{M+1} needed for S_{M+1}
    % S is (M+2) x (M+3)

    % one at top left
    sparseLowerDiagonal = speye(nCells+2, nCells+3);
    % one at bottom right
    sparseUpperDiagonal = sparse(1:nCells+2, 2:nCells+3,1, nCells+2, nCells+3);

    A = dog_math.getAMatrix(nBasisCpts);
    Phi1 = basis.getPhiMatrix(1, 1, nBasisCpts);
    Phi2 = basis.getPhiMatrix(-1, 1, nBasisCpts);

    C = -1/deltaX*Phi2;
    D = -1/deltaX*(A - Phi1);

    S = kron(sparseLowerDiagonal, C) + kron(sparseUpperDiagonal, D);
end
