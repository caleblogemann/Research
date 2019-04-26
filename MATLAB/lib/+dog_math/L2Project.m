function [q] = L2Project(q_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, time)
    q = zeros(num_cells, num_eqns, num_basis_cpts);
    dx = (b-a)/num_cells;
    x_i = @(i) a + (i - 1)*dx + dx/2;
    x = @(xi, i) x_i(i) + xi*dx/2;
    integrand_function = @(xi, i, k) q_function(x(xi, i), time)*basis.legendrePolynomial(k, xi);
    for k = 1:num_basis_cpts
        for i = 1:num_cells
            q(i, :, k) = 1/2*dog_math.gaussQuadrature(@(xi) integrand_function(xi, i, k), quad_order, -1, 1);
        end
    end
end
