function [R] = getRMatrix(nCells, nBasisCpts, deltaX, bc)
    % R_i = -1/deltaX*(A + Phi(-1, -1))Q_i + 1/deltaX*Phi(1, -1)Q_{i+1}
    % R_i = V*Q_i + W*Q_{i+1}
    % extrapolation
    % R = [ V+W 0 0   0 ]
    %     [ V+W 0 0   0 ]
    %     [   V W 0   0 ]
    %     [   0 V W   0 ]
    %     [   0 0 V   W ]
    %     [   0 0 0 V+W ]
    %     [   0 0 0 V+W ]
    %
    % periodic 
    % R = [ 0 0 V W ]
    %     [ W 0 0 V ]
    %     [ V W 0 0 ]
    %     [ 0 V W 0 ]
    %     [ 0 0 V W ]
    %     [ W 0 0 V ]
    %     [ V W 0 0 ]
    %
    %     Added first two rows to compute R_{-1} and R_0
    %     Added last row for R_{M+1}
    %     Don't add any columns so final matrix is M x M
    %     Use boundary conditions instead

    % one at (3,1)
    sparseLowerDiagonal = sparse(3:nCells+2, 1:nCells, 1, nCells+3, nCells);
    % one at (2, 1)
    sparseUpperDiagonal = sparse(2:nCells+1, 1:nCells, 1, nCells+3, nCells);

    A = dog_math.getAMatrix(nBasisCpts);
    Phi1 = basis.getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = basis.getPhiMatrix(1, -1, nBasisCpts);

    V = -1/deltaX*(A + Phi1);
    W = 1/deltaX*Phi2;

    R = kron(sparseLowerDiagonal, V) + kron(sparseUpperDiagonal, W);

    % Boundary Conditions
    nRows = (nCells+3)*nBasisCpts;
    nCols = nCells*nBasisCpts;
    firstRow = 1:nBasisCpts;
    firstCol = 1:nBasisCpts;
    secondRow = nBasisCpts+1:2*nBasisCpts;
    secondCol = nBasisCpts+1:2*nBasisCpts;
    secondLastRow = (nRows - 2*nBasisCpts + 1):(nRows - nBasisCpts);
    secondLastCol = (nCols - 2*nBasisCpts + 1):(nCols - nBasisCpts);
    lastRow = (nRows - nBasisCpts + 1):nRows;
    lastCol = (nCols - nBasisCpts + 1):nCols;
    if (strcmp(bc, 'extrapolation'))
        R(firstRow, firstCol) = V+W;
        R(secondRow, firstCol) = V+W;
        R(secondLastRow, lastCol) = V+W;
        R(lastRow, lastCol) = V+W;
    elseif (strcmp(bc,'periodic'))
        R(firstRow, secondLastCol) = V;
        R(firstRow, lastCol) = W;
        R(secondRow, lastCol) = V;
        R(secondLastRow, firstCol) = W;
        R(lastRow, firstCol) = V;
        R(lastRow, secondCol) = W;
    end
end
