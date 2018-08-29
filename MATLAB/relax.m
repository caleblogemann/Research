function [q] = relax(q, rhs, matrixFunction)
    global num_relaxations weight
    [num_elems, num_basis_cpts] = size(q);
    A = matrixFunction(num_elems);
    rhsVector = getVector(rhs);
    qVector = getVector(q);
    qVector = weighted_jacobi_relaxation(A, qVector, rhsVector, num_relaxations, weight );
    q = getMatrix(qVector, num_elems, num_basis_cpts);
%     D = diag(diag(A));
%     LPlusU = A - D;
% 
%     for i = 1:num_relaxations
%         qVector = (weight*(D\(rhsVector - LPlusU*qVector)) + (1 - weight)*qVector);
%     end
    
end
