function [A] = LDGThinFilmBE(num_elems, num_basis_cpts, original_q, deltaT, a, b)
    deltaX = (b - a)/num_elems;
    [original_num_elems, ~] = size(original_q);
    num_restrictions = log2(original_num_elems/num_elems);
    restricted_original_q = original_q;
    for i=1:num_restrictions
        restricted_original_q = restrict(restricted_original_q);
    end
    
    LDGMatrix = getLDGThinFilmMatrix(restricted_original_q, num_elems, num_basis_cpts, deltaX);
    
    A = eye(num_elems*num_basis_cpts) - deltaT*LDGMatrix;
end