function [LDG] = getLDGThinFilmMatrix(q, nCells, nBasisCpts, deltaX)
    R = getRMatrix(nCells, nBasisCpts, deltaX);
    S = getSMatrix(nCells, nBasisCpts, deltaX);
    U = getUMatrix(q, nCells, nBasisCpts, deltaX);
    Q = getQMatrix(nCells, nBasisCpts, deltaX);

    LDG = Q*U*S*R;
end
