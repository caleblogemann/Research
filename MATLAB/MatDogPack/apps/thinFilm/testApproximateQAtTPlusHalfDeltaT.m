exact_solution_function = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;
forcingFunction = @(x, t) 2.0*exp(-40.0*(t+30.0*(x-0.5).^2)).*...
    (648.0*exp(20.0*(t+30.0*(x-0.5).^2)).*...
    (14551 - 118200*x + 358200*x.^2 - 480000*x.^3 + 240000*x.^4)+...
    1296*exp(10.0*(t+30.0*(x-0.5).^2)).*...
    (21901 - 177600*x + 537600*x.^2 - 720000*x.^3 + 360000*x.^4)+...
    864.0*(29251 - 237000*x + 717000*x.^2 - 960000*x.^3 + 480000*x.^4)+...
    exp(30.0*(t+30.0*(x-0.5).^2)).*...
    (777707 - 6350400*x + 19310400*x.^2 - 25920000*x.^3 + 12960000*x.^4));

a = 0.0;
b = 1.0;

num_cells = 40;
num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;
dx = (b - a)/num_cells;

t = 0.0;

x = (a + dx/2):dx:(b-dx/2);
forcing_function = @(t) forcingFunction(x',t);

q_DG = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, 0.0);
% qOnes = dog_math.projectQ(@(x) ones(size(x)), nBasisCpts, nBasisCpts, num_cells, a, b);
qExact = @(t) dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, t);

for dt = (1/2).^(1:10)
    q_t_plus_half_delta_t = ApproximateQAtTPlusHalfDeltaT(q_DG, qExact, t, dt, dx, forcing_function);
    q_exact_t_plus_half_delta_t = qExact(t + 0.5*dt);
    error = sqrt(sum(sum((q_exact_t_plus_half_delta_t-q_t_plus_half_delta_t).^2)))/sqrt(sum(sum((q_exact_t_plus_half_delta_t).^2)))
end