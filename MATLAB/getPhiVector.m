function [phi] = getPhiVector(xi, nBasisCpts)
    % return column vector
    phi = sparse(arrayfun(@(n) legendrePolynomial(n, xi), 1:nBasisCpts)');
end
