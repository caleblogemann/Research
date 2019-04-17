function [new_solution] = IRK2(previous_solution, get_discretization_matrix, forcing_function, dt, old_time)
    [system_size, ~] = size(previous_solution);
    I = eye(system_size);

    A = get_discretization_matrix(old_time);

    % first stage
    rhs = (I + 0.25*dt*A)*previous_solution + 0.25*dt*(forcing_function(old_time) + forcing_function(old_time + 0.5*dt));
    qstar = (I - 0.25*dt*A)\rhs;

    % second stage
    rhs = (4*qstar - previous_solution + dt*forcing_function(old_time + dt));
    new_solution = (3*I - dt*A)\rhs;
end
