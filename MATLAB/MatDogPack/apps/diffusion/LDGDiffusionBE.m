function [A] = LDGDiffusionBE(num_elems, num_basis_cpts, deltaT, a, b, diffusivity, bc)
    deltaX = (b - a)/num_elems;

    LDGMatrix = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);

    A = eye(num_elems*num_basis_cpts) - deltaT*LDGMatrix;
end
