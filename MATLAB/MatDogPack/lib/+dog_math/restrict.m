function [result] = restrict(q)
    [num_elems, num_basis_cpts] = size(q);
    L = zeros(num_basis_cpts);
    R = zeros(num_basis_cpts);
    if (num_basis_cpts ==1)
        L(1,1)=0.5000000000000000;
        R(1,1)=0.5000000000000000;
    elseif (num_basis_cpts == 2)
        L(1,1)=0.5000000000000000;
        L(2,1)=-0.4330127018922193;
        L(2,2)=0.2500000000000000;
        R(1,1)=0.5000000000000000;
        R(2,1)=0.4330127018922193;
        R(2,2)=0.2500000000000000;
    else
        error('This number of basis components is not supported');
    end

    result = zeros(num_elems/2, num_basis_cpts);
    for i = 1:num_elems/2
        result(i,:) = (L*q(2*i-1,:)' + R*q(2*i,:)')';
    end
end
