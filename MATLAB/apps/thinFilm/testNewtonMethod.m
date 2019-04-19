% Test Newton iteration
% exact_solution_function = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;
% forcingFunction = @(x, t) 2.0*exp(-40.0*(t+30.0*(x-0.5).^2)).*...
%     (648.0*exp(20.0*(t+30.0*(x-0.5).^2)).*...
%     (14551 - 118200*x + 358200*x.^2 - 480000*x.^3 + 240000*x.^4)+...
%     1296*exp(10.0*(t+30.0*(x-0.5).^2)).*...
%     (21901 - 177600*x + 537600*x.^2 - 720000*x.^3 + 360000*x.^4)+...
%     864.0*(29251 - 237000*x + 717000*x.^2 - 960000*x.^3 + 480000*x.^4)+...
%     exp(30.0*(t+30.0*(x-0.5).^2)).*...
%     (777707 - 6350400*x + 19310400*x.^2 - 25920000*x.^3 + 12960000*x.^4));

exact_solution_function = @(x, t) 0.1*sin(2*pi*(x - t)) + 0.15;
forcingFunction = @(x,t) -0.628319*cos(2*pi*(x-t)) ... 
    - 46.7564*cos(2*pi*(x-t)).^2.*(0.15 + 0.1*sin(2*pi*(x-t))).^2 ...
    + 155.855*(0.15 + 0.1*sin(2*pi*(x-t))).^3.*sin(2*pi*(x-t));

a = 0.0;
b = 1.0;

num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

time_order = 1;
initial_time = 0.0;
final_time = 0.1;

num_doublings = 3;
err = zeros(1, num_doublings);
initial_num_cells = 100;
for i = 1:num_doublings
    
    num_cells = initial_num_cells*2^(i-1);
    deltaX = (b - a)/num_cells;

    x = (a + 0.5*deltaX):deltaX:(b-0.5*deltaX);
    forcing_function = @(t) forcingFunction (x', t);
 
    deltaT = deltaX;
    num_time_steps = final_time/deltaT;
    
    q_FD = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, initial_time);
    q_FD_exact = @(t) dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, t);
    I = eye(num_cells);
    
    residual_array = zeros(num_time_steps,1);
    iteration_array = zeros(num_time_steps, 1);
    error_array = zeros(num_time_steps, 1);
    for n = 1:num_time_steps
        old_time = initial_time + (n-1)*deltaT;
        q_FD_old = q_FD;
        F = @(q) q - deltaT*getFDThinFilmMatrix(q, deltaX)*q - q_FD_old - deltaT*forcing_function(old_time+deltaT);
        J = @(q) I - deltaT*getFDThinFilmJacobian(q, deltaX);
        res = @(q) norm(F(q))/norm(q);

        % initial guess
        q_FD = q_FD_old;
        
        residual = res(q_FD);
        iter = 0;
        max_num_newton_iterations = 10;
        while(residual > 1e-5 && iter <= max_num_newton_iterations)
            % x_{n+1} = x_n - J(x_n)^{-1} F(x_n)

            delta = J(q_FD)\F(q_FD);
            q_FD = q_FD - delta;%J(q_FD)\F(q_FD);
            residual = res(q_FD);
            iter = iter + 1;
        end
        residual_array(n) = residual;
        error_array(n) = dog_math.ComputeError(q_FD, exact_solution_function, a, b, old_time+deltaT);
        iteration_array(n) = iter;
    end
    err(i) = dog_math.ComputeError(q_FD, exact_solution_function, a, b, final_time);
    plot(x, q_FD, x, exact_solution_function(x, final_time));
    pause();
    plot(residual_array);
    pause();
    plot(error_array);
    pause();
    plot(iteration_array);
    pause();
end
deltaXArray = (b-a)./(initial_num_cells*2.^(0:num_doublings-1));
log(err(1:end-1)./err(2:end))./log(deltaXArray(1:end-1)./deltaXArray(2:end))