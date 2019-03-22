function [q] = projectQ(qFunction, quadOrder, nBasisCpts, nCells, a, b)
    q = zeros(nCells, nBasisCpts);
    dx = (b-a)/nCells;
    x_i = @(i) a + (i - 1)*dx + dx/2;
    x = @(xi, i) x_i(i) + xi*dx/2;
    for k = 1:nBasisCpts
        for i = 1:nCells
            q(i, k) = 1/2*gaussQuadrature(quadOrder, @(xi) qFunction(x(xi, i))*legendrePolynomial(k, xi));
        end
    end
end
