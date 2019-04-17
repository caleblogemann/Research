function [phi] = getPhiVector(xi, nBasisCpts)
    % return column vector
    phi = sparse(arrayfun(@(n) basis.legendrePolynomial(n, xi), 1:nBasisCpts)');
end
