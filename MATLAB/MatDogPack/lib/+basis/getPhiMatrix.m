function [Phi] = getPhiMatrix(xi1, xi2, nBasisCpts)
    Phi = basis.getPhiVector(xi1, nBasisCpts)*basis.getPhiVector(xi2, nBasisCpts)';
end
