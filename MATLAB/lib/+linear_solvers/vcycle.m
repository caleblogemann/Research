function [q] = vcycle(q, rhs, matrixFunction)
    global vcycle_depth
    
    [num_elems, num_basis_cpts] = size(q);

    q_array = cell(vcycle_depth+1,1);
    rhs_array = cell(vcycle_depth+1,1);
    q_array(1) = {q};
    rhs_array(1) = {rhs};

    temp_n_elems = num_elems/2;
    for i = 2:vcycle_depth+1
        q_array(i) = {zeros(temp_n_elems, num_basis_cpts)};
        rhs_array(i) = {zeros(temp_n_elems, num_basis_cpts)};
        temp_n_elems = temp_n_elems/2;
    end

    current_num_elems = num_elems;
    for i = 1:vcycle_depth
        q_array(i) = {relax(q_array{i}, rhs_array{i}, matrixFunction)};
        Aq = matrixFunction(current_num_elems)*getQVector(q_array{i});
        residual = getQVector(rhs_array{i}) - Aq;
        rhs_array{i+1} = restrict(getMatrix(residual, current_num_elems, num_basis_cpts));
        current_num_elems = current_num_elems/2;
    end

    % solve exactly
    A = matrixFunction(current_num_elems);
    q_array{vcycle_depth+1} = getMatrix(A\getVector(rhs_array{vcycle_depth+1}), current_num_elems, num_basis_cpts);

    for i = vcycle_depth:-1:1
        current_num_elems = 2*current_num_elems;
        err = interpolate(q_array{i+1});
        q_array{i} = q_array{i} + err;
        q_array(i) = {relax(q_array{i}, rhs_array{i}, matrixFunction)};
    end

    q = q_array{1};
end
