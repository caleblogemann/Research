% Test forcing function
exact_solution_function = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;

forcingFunction = @(x, t) 2.0*exp(-40.0*(t+30.0*(x-0.5).^2)).*...
    (648.0*exp(20.0*(t+30.0*(x-0.5).^2)).*...
    (14551 - 118200*x + 358200*x.^2 - 480000*x.^3 + 240000*x.^4)+...
    1296*exp(10.0*(t+30.0*(x-0.5).^2)).*...
    (21901 - 177600*x + 537600*x.^2 - 720000*x.^3 + 360000*x.^4)+...
    864.0*(29251 - 237000*x + 717000*x.^2 - 960000*x.^3 + 480000*x.^4)+...
    exp(30.0*(t+30.0*(x-0.5).^2)).*...
    (777707 - 6350400*x + 19310400*x.^2 - 25920000*x.^3 + 12960000*x.^4));

a = 0;
b = 1;
num_cells = 400;
num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

deltaX = (b - a)/num_cells;
x = (a+deltaX/2):deltaX:(b - deltaX/2);

old_time = 0.0;
deltaT = 0.00001;

q_DG_n = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, old_time);
q_FD_n = dog_math.getQFD(q_DG_n);
q_DG_np1 = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, old_time+deltaT);
q_FD_np1 = dog_math.getQFD(q_DG_np1);

F = FDThinFilmOperator(q_FD_np1, deltaX);

q_t = (q_FD_np1 - q_FD_n)/deltaT;
disp(norm(q_t - F - forcingFunction(x, old_time + deltaT))/norm(q_FD_np1));