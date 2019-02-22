function [q, num_iterations] = multigrid(q, rhs, matrixFunction)
    global tolerance
    global max_iterations
    num_iterations = 0;
    [num_elems, ~] = size(q);
    A = matrixFunction(num_elems);
    qVector = getVector(q);
    Aq = A*qVector;
    residual = Aq - getQVector(rhs);
    rr = norm(residual, 'fro');
    while (num_iterations < max_iterations && rr > tolerance)
        num_iterations = num_iterations + 1;
        q = vcycle(q, rhs, matrixFunction);
        Aq = A*getQVector(q);
        residual = Aq - getQVector(rhs);
        rr = norm(residual, 'fro');
    end
    if (num_iterations == max_iterations)
        disp('MultiGrid Did not converge');
        disp('rr = ');
        disp(rr);
    end
end
