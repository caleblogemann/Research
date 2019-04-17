function [ new_solution ] = TimeStep( previous_solution, get_discretization_matrix, forcing_function, time_order, dt, old_time )
    % Solve q_t = A q + f
    % get_discretization_matrix gives A and should be a function of t
    % forcing function gives vector f and should also be function of t
    % dt is time step
    % old_time is time at previous_solution
    if(time_order == 1)
        new_solution = time_stepper.ImplicitRungeKutta.BackwardEuler(previous_solution, get_discretization_matrix, forcing_function, dt, old_time);
    elseif (time_order == 2)
        new_solution = time_stepper.ImplicitRungeKutta.IRK2(previous_solution, get_discretization_matrix, forcing_function, dt, old_time);
    end
end

