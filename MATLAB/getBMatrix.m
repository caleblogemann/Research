function [B] = getBMatrix(nBasisCpts, q, i)
    % B_{kl} = int{-1}{1}{q^3_i phi^k(xi) phi^l_{xi}(xi)}{xi}
    B = sparse(nBasisCpts, nBasisCpts);
    qCubedFunction = @(xi) sum(arrayfun(@(m) q(i, m)*legendrePolynomial(m, xi), 1:nBasisCpts))^3;
    for k = 1:nBasisCpts
        for l = 2:nBasisCpts
            integrandFunction = @(xi) qCubedFunction(xi)*legendrePolynomial(k, xi)*legendrePolynomialDerivative(l, xi);
            B(k, l) = gaussQuadrature(nBasisCpts+1, integrandFunction);
        end
    end
end
