exact_solution_function = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;

a = 0.0;
b = 1.0;

num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

err = [];
n = 10;
for dx = (0.5).^(1:n)
    num_cells = (b - a)/dx;
    q_DG = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, 0.0);
    q_FD = dog_math.getQFD(q_DG);
    q_DG = dog_math.getQDG(q_FD);
    err = [err, dog_math.ComputeError(q_DG, exact_solution_function, a, b, 0.0)];
end

log(err(1:end-1)./err(2:end))./log(((0.5).^(1:n-1))./((0.5).^(2:n)))