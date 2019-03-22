function [result] = interpolate(q)
    [num_elems, num_basis_cpts] = size(q);
    L = zeros(num_basis_cpts);
    R = zeros(num_basis_cpts);
    if (num_basis_cpts ==1)
        L(1,1)=1.000000000000000;
        R(1,1)=1.000000000000000;
    elseif (num_basis_cpts == 2)
        L(1,1)=1.000000000000000;
        L(1,2)=-0.8660254037844386;
        L(2,2)=0.5000000000000000;
        R(1,1)=1.000000000000000;
        R(1,2)=0.8660254037844386;
        R(2,2)=0.5000000000000000;
    else
        error('This number of basis components is not supported');
    end
    result = zeros(num_elems*2, num_basis_cpts);
    for i = 1:num_elems
        result(2*i-1,:) = (L*q(i,:)')';
        result(2*i,:) = (R*q(i,:)')';
    end
end
