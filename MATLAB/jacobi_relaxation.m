function [ sol ] = jacobi_relaxation( A, sol, rhs, num_relaxations )
    % q_{n+1} = D^{-1}(rhs - (L+U)sol_n)
    D = diag(diag(A));
    LPlusU = A - D;
    for i = 1:num_relaxations
        sol = D\(rhs - LPlusU*sol);
    end
end

