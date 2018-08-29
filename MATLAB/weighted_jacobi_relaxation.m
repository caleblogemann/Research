function [ sol ] = weighted_jacobi_relaxation( A, sol, rhs, num_relaxations, weight )
    % q_{n+1} = weight*D^{-1}(rhs - (L+U)sol_n) + (1 - weight)*sol_n
    D = diag(diag(A));
    LPlusU = A - D;
    for i = 1:num_relaxations
        sol = weight*(D\(rhs - LPlusU*sol)) + (1 - weight)*sol;
    end
end