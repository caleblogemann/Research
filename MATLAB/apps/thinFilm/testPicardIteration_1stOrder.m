% Test picard iteration
% exact_solution_function = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;
% forcingFunction = @(x, t) 2.0*exp(-40.0*(t+30.0*(x-0.5).^2)).*...
%     (648.0*exp(20.0*(t+30.0*(x-0.5).^2)).*...
%     (14551 - 118200*x + 358200*x.^2 - 480000*x.^3 + 240000*x.^4)+...
%     1296*exp(10.0*(t+30.0*(x-0.5).^2)).*...
%     (21901 - 177600*x + 537600*x.^2 - 720000*x.^3 + 360000*x.^4)+...
%     864.0*(29251 - 237000*x + 717000*x.^2 - 960000*x.^3 + 480000*x.^4)+...
%     exp(30.0*(t+30.0*(x-0.5).^2)).*...
%     (777707 - 6350400*x + 19310400*x.^2 - 25920000*x.^3 + 12960000*x.^4));

amplitude = 0.2;
steepness = -12.0;
wavespeed = 1.0;
center = 5.0;
offset = 0.1;
final_time = 1.0;
a = 0;
b = 10.0;
exact_solution_function = @(x, t) amplitude*exp(steepness*(x - wavespeed*t - center).^2) + offset;
q_x = @(x, t) amplitude*2.0*steepness*(x - wavespeed*t - center)...
    .*exp(steepness*(x - wavespeed*t - center).^2.0);
q_t = @(x, t) -1.0*wavespeed*q_x(x, t);
q_xxx = @(x, t) amplitude*(12.0*exp(steepness*(x - wavespeed*t - center).^2.0)...
    .*steepness^2.0.*(x - wavespeed*t - center) ...
    + 8.0*exp(steepness*(x - wavespeed*t - center).^2.0)...
    .*steepness^3.0.*(x - wavespeed*t - center).^3.0);
q_xxxx = @(x, t) amplitude*(12.0*exp(steepness*(x - wavespeed*t - center).^2.0).*steepness.^2.0 ...
    + 48.0*exp(steepness*(x - wavespeed*t - center).^2.0).*steepness^3.0.*(x - wavespeed*t - center).^2.0 ...
    + 16.0*exp(steepness*(x - wavespeed*t - center).^2.0).*(steepness*(x - wavespeed*t -center)).^4.0);

% amplitude = 0.1;
% wavespeed = 1.0;
% wavenumber = 1.0;
% offset = 0.15;
% final_time = 1.0;
% a = 0;
% b = 1.0;
% exact_solution_function = @(x, t) amplitude*sin(2*pi*wavenumber*(x - wavespeed*t)) + offset;
% q_t = @(x, t) -1.0*amplitude*2*pi*wavenumber*wavespeed*cos(2*pi*wavenumber*(x - wavespeed*t));
% q_x = @(x, t) amplitude*2*pi*wavenumber*cos(2*pi*wavenumber*(x - wavespeed*t));
% q_xxx = @(x, t) -amplitude*(2*pi*wavenumber)^3*cos(2*pi*wavenumber*(x - wavespeed*t));
% q_xxxx = @(x, t) amplitude*(2*pi*wavenumber)^4*sin(2*pi*wavenumber*(x - wavespeed*t));

% exact_solution_function = @(x, t) (0.3)*(x - t <= 0.5) + (0.1)*(x - t > 0.5);
% forcingFunction = @(x, t) zeros(size(x));

forcingFunction = @(x,t) q_t(x, t) + 3.0*exact_solution_function(x, t).^2.0.*q_x(x, t).*q_xxx(x, t) + exact_solution_function(x, t).^3.0.*q_xxxx(x, t);

forcing = true;
pausing = true;
cfl = 1.0;

max_num_iterations = 2;
tolerance = 1e-5;

num_eqns = 1;
num_basis_cpts = 2;
quad_order = 2;

time_order = 1;
initial_time = 0.0;

num_doublings = 4;
err = zeros(1, num_doublings);
initial_num_cells = 200;
for i = 1:num_doublings
    
    num_cells = initial_num_cells*2^(i-1);
    deltaX = (b - a)/num_cells;

    x = (a + 0.5*deltaX):deltaX:(b-0.5*deltaX);
    if(forcing)
        forcing_function = @(t) forcingFunction (x', t);
    else
        forcing_function = @(t) zeros(size(x'));
    end
 
    deltaT = cfl*deltaX;
    num_time_steps = round(final_time/deltaT);
    deltaT = final_time/num_time_steps;
    
    q_FD = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, initial_time);
    q_FD_exact = @(t) dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, 1, a, b, t);
    I = eye(num_cells);
    
    residual_array_1 = zeros(num_time_steps,1);
    iteration_array_1 = zeros(num_time_steps, 1);
    error_array = zeros(num_time_steps, 1);
    for n = 1:num_time_steps
        old_time = initial_time + (n-1)*deltaT;
        q_FD_old = q_FD;
        
        % 1st stage
        F = @(q) q - deltaT*FDThinFilmOperator(q, deltaX) - q_FD_old - deltaT*forcing_function(old_time+deltaT);
        res = @(q) norm(F(q))/norm(q);
  
        residual = res(q_FD);
        iter = 0;
        
        rhs = q_FD_old + deltaT*forcing_function(old_time+deltaT);
        % initial guess
        q_FD = q_FD_old;
        while(residual > tolerance && iter < max_num_iterations)
                A = getFDThinFilmMatrix(q_FD, deltaX);
                q_FD = (I - deltaT*A)\rhs;
                %q_FD(q_FD <= 0) = 1e-5;
                residual = res(q_FD);
                iter = iter+1;
        end
        residual_array_1(n) = residual;
        iteration_array_1(n) = iter;
        error_array(n) = dog_math.ComputeError(dog_math.getQDG(q_FD), exact_solution_function, a, b, old_time+deltaT);
    end
    err(i) = dog_math.ComputeError(dog_math.getQDG(q_FD), exact_solution_function, a, b, final_time);
    plot(x, q_FD, x, exact_solution_function(x, final_time));
    title('Approximate Solution');
    xlabel('x');
    if(pausing)
        pause();
    end
    plot(1:num_time_steps, residual_array_1);
    title('Residuals');
    xlabel('Time Steps');
    if(pausing)
        pause();
    end
    plot(error_array);
    title('Errors');
    xlabel('Time Steps');
    if(pausing)
        pause();
    end
    plot(1:num_time_steps, iteration_array_1);
    title('Number of Iterations');
    xlabel('Time Steps');
    if(pausing)
        pause();
    end
end
deltaXArray = (b-a)./(initial_num_cells*2.^(0:num_doublings-1));
log(err(1:end-1)./err(2:end))./log(deltaXArray(1:end-1)./deltaXArray(2:end))