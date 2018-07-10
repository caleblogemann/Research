function [S] = getSMatrix(nCells, nBasisCpts, deltaX)
    % S_i = -1/deltaX*(A - Phi(1, 1))R_i - 1/deltaX*Phi(-1, 1)R_{i-1}
    I = speye(nCells, nCells);
    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(1, 1, nBasisCpts);
    Phi2 = getPhiMatrix(-1, 1, nBasisCpts);

    B = -1/deltaX*(A - Phi1);
    C = -1/deltaX*Phi2;

    sparseLowerDiagonal = sparse(2:nCells, 1:nCells-1, 1, nCells, nCells);

    S = kron(I, B) + kron(sparseLowerDiagonal, C);
    % ghost cell at left endpoint
    S(1:nBasisCpts, 1:nBasisCpts) = S(1:nBasisCpts, 1:nBasisCpts) + C;
end
