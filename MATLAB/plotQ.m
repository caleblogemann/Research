function [] = plotQ(q, a, b)
    [nCells, nBasisCpts] = size(q);
    dx = (b - a)/nCells;
    i = @(x) floor((x - a)/dx)+1;
    x_i = @(i) a + (i - 1)*dx + dx/2;
    xi = @(x, i) 2/dx*(x - x_i(i));
    qFunction = @(x) arrayfun(@(xj) sum(arrayfun(@(l) q(i(xj), l)*legendrePolynomial(l, xi(xj, i(xj))),1:nBasisCpts)) , x);
    x = linspace(a, b-eps, nCells*10);
    plot(x, qFunction(x), 'k');
end

