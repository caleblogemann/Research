function [A] = LDGHyperDiffusionBE(num_elems, num_basis_cpts, deltaT, a, b, diffusivity, bc)
    deltaX = (b - a)/num_elems;
    LDGMatrix = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
    A = eye(num_elems) - deltaT*LDGMatrix;
end
