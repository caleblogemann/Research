function [R] = getRMatrixDiffusion(nCells, nBasisCpts, deltaX, bc)
    % R_i = -1/deltaX*(A + Phi(-1, -1))Q_i + 1/deltaX*Phi(1, -1)Q_{i+1}
    % R_i = V*Q_i + W*Q_{i+1}
    % extrapolation
    % R = [ V+W 0 0   0 ]
    %     [   V W 0   0 ]
    %     [   0 V W   0 ]
    %     [   0 0 V   W ]
    %     [   0 0 0 V+W ]
    % periodic 
    % R = [ W 0 0 V ]
    %     [ V W 0 0 ]
    %     [ 0 V W 0 ]
    %     [ 0 0 V W ]
    %     [ W 0 0 V ]
    %
    %     Added first row to compute R_0
    %     Don't add any columns so final matrix is M x M
    %     Use boundary conditions instead

    % one at (2,1)
    sparseLowerDiagonal = sparse(2:nCells+1, 1:nCells, 1, nCells+1, nCells);
    % one at (1, 1)
    sparseUpperDiagonal = sparse(1:nCells, 1:nCells, 1, nCells+1, nCells);

    A = getAMatrix(nBasisCpts);
    Phi1 = getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = getPhiMatrix(1, -1, nBasisCpts);

    V = -1/deltaX*(A + Phi1);
    W = 1/deltaX*Phi2;

    R = kron(sparseLowerDiagonal, V) + kron(sparseUpperDiagonal, W);

    % Boundary Conditions
    nRows = (nCells+1)*nBasisCpts;
    nCols = nCells*nBasisCpts;
    firstRow = 1:nBasisCpts;
    firstCol = 1:nBasisCpts;
    lastRow = (nRows - nBasisCpts + 1):nRows;
    lastCol = (nCols - nBasisCpts + 1):nCols;
    if (strcmp(bc, 'extrapolation'))
        R(firstRow, firstCol) = V+W;
        R(lastRow, lastCol) = V+W;
    elseif (strcmp(bc, 'periodic'))
        R(firstRow, lastCol) = V;
        R(lastRow, firstCol) = W;
    end
end
