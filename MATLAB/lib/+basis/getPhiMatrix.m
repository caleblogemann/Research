function [Phi] = getPhiMatrix(xi1, xi2, nBasisCpts)
    Phi = getPhiVector(xi1, nBasisCpts)*getPhiVector(xi2, nBasisCpts)';
end
