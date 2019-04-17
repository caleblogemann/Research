function [err] = ComputeError(approximate_solution, exact_solution_function, a, b, time)
    % approximate solution in DG form
    % exact_solution_function - function of (x, t)

    % asumming in 1 dimension

    [num_cells, num_eqns, num_basis_cpts] = size(approximate_solution);
    quad_order = 10;
    exact_solution = dog_math.L2Project(exact_solution_function, quad_order, num_cells, num_eqns, num_basis_cpts + 1, a, b, time);

    difference = exact_solution;
    difference(:, :, 1:num_basis_cpts) = exact_solution(:, :, 1:num_basis_cpts) - approximate_solution;
    
    exact_solution_norm = sqrt(sum(sum(exact_solution.^2)));
    difference_norm = sqrt(sum(sum(difference.^2)));
    err = difference_norm/exact_solution_norm;
end
