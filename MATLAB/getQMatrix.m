function [Q] = getQMatrix(nCells, nBasisCpts, deltaX)
    % Q_i = 1/deltaX*(A - Phi(1, 1))U_i + 1/deltaX*Phi(-1, 1)U_{i-1}
    I = speye(nCells, nCells);
    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(1, 1, nBasisCpts);
    Phi2 = getPhiMatrix(-1, 1, nBasisCpts);

    C = 1/deltaX*(A - Phi1);
    D = 1/deltaX*Phi2;

    sparseLowerDiagonal = sparse(2:nCells, 1:nCells-1, 1, nCells, nCells);

    Q = kron(I, C) + kron(sparseLowerDiagonal, D);
    % ghost cell at left endpoint
    Q(1:nBasisCpts, 1:nBasisCpts) = Q(1:nBasisCpts, 1:nBasisCpts) + D;
end
