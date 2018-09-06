function [U] = getUMatrix(q, nCells, nBasisCpts, deltaX)
    % U_i = 1/deltaX*(B_i - q^-_{i+1/2}^3*Phi(1, 1) + (q^+_{i-1/2}^3 - \hat{q}^3_{i-1/2})*Phi(-1,-1))S_i
    % + 1/deltaX*\hat{q}^3_{i+1/2}*Phi(1, -1)S_{i+1}
    Phi0 = getPhiMatrix(1, 1, nBasisCpts);
    Phi1 = getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = getPhiMatrix(1, -1, nBasisCpts);

    % left side at xi = 1
    phiVectorMinus = arrayfun(@(l) legendrePolynomial(l,  1), 1:nBasisCpts);
    % right side at xi = -1
    phiVectorPlus  = arrayfun(@(l) legendrePolynomial(l, -1), 1:nBasisCpts);
    % qminus(i) = q^-_{i-1/2}^3
    qminus = ([q(1,:); q]*phiVectorMinus').^3;
    % qplus(i) = q^+_{i-1/2}^3
    qplus  = ([q; q(end,:)]*phiVectorPlus').^3;
    % qhat(i) = (q^-_{i-1/2}^3 + q^+_{i-1/2}^3)/2
    qhat   = (qminus+qplus)/2;
    C = @(i) 1/deltaX*(getBMatrix(nBasisCpts, q, i) - qminus(i+1)*Phi0 + (qplus(i) - qhat(i))*Phi1);
    D = @(i) 1/deltaX*qhat(i+1)*Phi2;

    U = sparse(nCells*nBasisCpts, nCells*nBasisCpts);
    for i = 1:nCells
        if(i < nCells)
            U((nBasisCpts*(i-1)+1):nBasisCpts*i,(nBasisCpts*(i-1)+1):nBasisCpts*i) = C(i);
            U((nBasisCpts*(i-1)+1):nBasisCpts*i,(nBasisCpts*(i)+1):nBasisCpts*(i+1)) = D(i);
        else
            % ghost cell at right endpoint
            U((nBasisCpts*(i-1)+1):nBasisCpts*i,(nBasisCpts*(i-1)+1):nBasisCpts*i) = C(i) + D(i);
        end
    end
end
