function [new_solution] = BackwardEuler(previous_solution, get_discretization_matrix, forcing_function, deltaT, old_time)
    [system_size, ~] = size(previous_solution);
    I = eye(system_size);
    
    rhs = previous_solution + deltaT*forcing_function(old_time + deltaT);
    
    A = get_discretization_matrix(old_time);
    new_solution = (I - deltaT*A)\rhs;
end
