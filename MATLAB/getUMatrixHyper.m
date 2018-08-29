function [U] = getUMatrixHyper(nCells, nBasisCpts, deltaX)
    % U_i = -1/deltaX*(A + Phi(-1, -1))S_i + 1/deltaX*Phi(1, -1)S_{i+1}
    I = speye(nCells, nCells);
    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = getPhiMatrix(1, -1, nBasisCpts);

    B = -1/deltaX*(A + Phi1);
    C = 1/deltaX*Phi2;

    sparseUpperDiagonal = sparse(1:nCells-1, 2:nCells, 1, nCells, nCells);

    U = kron(I, B) + kron(sparseUpperDiagonal, C);
    % ghost cell at right endpoint
    U((end-nBasisCpts+1):end, (end-nBasisCpts+1):end) = U((end-nBasisCpts+1):end, (end-nBasisCpts+1):end) + C;
end
