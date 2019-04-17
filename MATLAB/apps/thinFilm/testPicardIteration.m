% Test picard iteration
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

num_cells = 100;
num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

deltaX = (b - a)/num_cells;

x = (a + 0.5*deltaX):deltaX:(b-0.5*deltaX);
forcing_function = @(t) forcingFunction (x', t);

old_time = 0.0;
deltaT = 0.1*deltaX;
time_order = 1.0;

q_DG = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts, a, b, old_time);
q_FD_old = dog_math.getQFD(q_DG);

F = @(q_FD) q_FD - deltaT*FDThinFilmOperator(q_FD, deltaX) - q_FD_old - deltaT*forcing_function(old_time+deltaT);
res = @(q_FD) norm(F(q_FD))/norm(q_FD);

err = [];
q_FD = q_FD_old;
n = 10;
for i = 1:n
    get_discretization_matrix = @(t) getFDThinFilmMatrix(q_FD, deltaX);
    % take a step with FD solution
    q_FD = time_stepper.ImplicitRungeKutta.TimeStep(q_FD_old, get_discretization_matrix, forcing_function, time_order, deltaT, old_time);
    % switch back to DG solution
    % q_DG = dog_math.getQDG(q_FD);
    disp(res(q_FD));
    %disp(dog_math.ComputeError(q_DG, exact_solution_function, a, b, old_time+dt));
end