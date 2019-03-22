function [ matrix ] = getMatrix( vector, num_elems, num_basis_cpts )
    matrix = zeros(num_elems, num_basis_cpts);
    for i = 1:num_elems
        for k = 1:num_basis_cpts
            matrix(i,k) = vector(num_basis_cpts*(i-1) + k);
        end
    end
end

