function [LDG] = getLDGThinFilmMatrix(q, nCells, nBasisCpts, deltaX, diffusivity, bc)
    R = getRMatrix(nCells, nBasisCpts, deltaX, bc);
    S = getSMatrix(nCells, nBasisCpts, deltaX);
    U = getUMatrix(q, nCells, nBasisCpts, deltaX);
    Q = getQMatrix(nCells, nBasisCpts, deltaX, diffusivity);

    LDG = Q*U*S*R;
end
