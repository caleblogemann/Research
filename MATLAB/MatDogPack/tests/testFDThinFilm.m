% for f = q
% fourth order operator
%lambda = 2*pi;
%exactSolution = @(x, t) exp(-1.0*lambda^4*t).*sin(lambda*x);
%exactSolution = @(x, t) 0.2*exp(-100.0*(x - 0.5).^2) + 0.1;
%forcingFunction = @(x) 8.0*exp(-100.0*(1.0 - 2.0*x).^2).*...
%    (2.0 + exp(25.0*(1.0 - 2.0*x).^2)).^2.*...
%    (exp(25.0*(1.0-2.0*x).^2).*(2203 - 18800*x + 58800*x.^2 - 80000*x.^3 + ...
%    40000*x.^4) + 2*(9253 - 77000*x + 237000*x.^2 - 320000*x.^3 + 160000*x.^4));

% time dependent exact solution
% exactSolution = @(x, t) 0.2*exp(-10*t).*exp(-100.0*(x - 0.5).^2) + 0.1;
% forcingFunction = @(x, t) exp(-40.0*(t+10.0*(x-0.5).^2)).*...
%     (288.0*exp(10.0*(t+10.0*(x-0.5).^2)).*...
%     (2301 - 19200*x + 59200*x.^2 - 80000*x.^3 + 40000*x.^4)+...
%     48.0*exp(20.0*(t+10.0*(x-0.5).^2)).*...
%     (4553 - 38200*x + 118200*x.^2 - 160000*x.^3 + 80000*x.^4)+...
%     2.0*exp(30.0*(t+10.0*(x-0.5).^2)).*...
%     (8811 - 75200*x + 235200*x.^2 - 320000*x.^3 + 160000*x.^4)+...
%     64*(9253 - 77000*x + 237000*x.^2 - 320000*x.^3 + 160000*x.^4));

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

num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

diffusivity = 1.0;
bc = 'periodic';

time_order = 1;
final_time = 0.1;

num_doublings = 3;
err = zeros(1, num_doublings);
initial_num_cells = 100;
for i = 1:num_doublings
    num_cells = initial_num_cells*2^(i-1);
    deltaX = (b-a)/num_cells;
    
    deltaT = 0.1*deltaX;
    nTimeSteps = final_time/deltaT;
    %final_time = nTimeSteps*deltaT;
    
    x = (a + deltaX/2):deltaX:(b-deltaX/2);
    forcing_function = @(t) forcingFunction(x',t);
    %forcing_function = @(t) zeros(num_cells, 1);

    q_FD = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, 0.0);
    q_FD_exact = @(t) dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, t);
    
    residual_array = zeros(nTimeSteps,1);
    iteration_array = zeros(nTimeSteps, 1);
    for n = 1:nTimeSteps
        old_time = (n-1)*deltaT;
        
        q_FD_old = q_FD;
        F = @(q_FD) q_FD - deltaT*FDThinFilmOperator(q_FD, deltaX) - q_FD_old - deltaT*forcing_function(old_time+deltaT);
        res = @(q_FD) norm(F(q_FD))/norm(q_FD);
        
        % initial guess
        q_FD = q_FD_old;
        residual = res(q_FD);
        iter = 0;
        % picard iteration, take time step=
        while (residual > .1 && iter < 100)
            get_discretization_matrix = @(t) getFDThinFilmMatrix(q_FD, deltaX);
            % take a step with FD solution
            q_FD = time_stepper.ImplicitRungeKutta.TimeStep(q_FD_old, get_discretization_matrix, forcing_function, time_order, deltaT, old_time);
            residual = res(q_FD);
            iter = iter+1;
        end
        iteration_array(n) = iter;
        residual_array(n) = residual;
    end

    err(i) = dog_math.ComputeError(q_FD, exact_solution_function, a, b, final_time);
    plot(x, q_FD, x, exact_solution_function(x, final_time));
    pause();
    plot(1:nTimeSteps, residual_array, 1:nTimeSteps, iteration_array);
    disp(i);
    pause()
end
deltaXArray = (b-a)./(initial_num_cells*2.^(0:num_doublings-1));
log(err(1:end-1)./err(2:end))./log(deltaXArray(1:end-1)./deltaXArray(2:end))
%log(error(1:end-1)./error(2:end))./log(deltaXArray(1:end-1)./deltaXArray(2:end))
