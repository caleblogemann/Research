function [ sol ] = gauss_seidel_relaxation( A, sol, rhs, num_relaxations )
    % q_{n+1} = weight*D^{-1}(rhs - (L+U)sol_n) + (1 - weight)*sol_n
    L = tril(A);
    U = A - L;
    for i = 1:num_relaxations
        sol = L\(rhs - U*sol);
    end
end

