function [A] = LDGHyperDiffusionBE(num_elems, num_basis_cpts, deltaT, a, b)
    deltaX = (b - a)/num_elems;

    LDGMatrix = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX);

    A = eye(num_elems) - deltaT*LDGMatrix;
end
