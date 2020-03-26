function [U] = getUMatrix(q, nCells, nBasisCpts, deltaX, bc)
    % U_i = 1/deltaX*(B_i - q^-_{i+1/2}^3*Phi(1, 1) + (q^+_{i-1/2}^3 - \hat{q}^3_{i-1/2})*Phi(-1,-1))S_i
    % + 1/deltaX*\hat{q}^3_{i+1/2}*Phi(1, -1)S_{i+1}
    % U_i = E_i S_i + F_i S_{i+1}
    % U = [ E F 0 0 0 0 ]
    %     [ 0 E F 0 0 0 ]
    %     [ 0 0 E F 0 0 ]
    %     [ 0 0 0 E F 0 ]
    %     [ 0 0 0 0 E F ]
    %     Add first row to compute U_0 for Q
    %     Add last column to use S_{M+1}
    %     Added first column to use S_{0}
    %     U is (M+1) x (M+2) because added row and 2 columns

    Phi0 = basis.getPhiMatrix(1, 1, nBasisCpts);
    Phi1 = basis.getPhiMatrix(-1, -1, nBasisCpts);
    Phi2 = basis.getPhiMatrix(1, -1, nBasisCpts);

    % left side at xi = 1
    phiVectorMinus = arrayfun(@(l) basis.legendrePolynomial(l,  1), 1:nBasisCpts);
    % right side at xi = -1
    phiVectorPlus  = arrayfun(@(l) basis.legendrePolynomial(l, -1), 1:nBasisCpts);
    switch bc
        case 'extrapolation'
            qExtended = [q(1,:,:); q(1,:,:); q; q(end,:,:)];
        case 'periodic'
            qExtended = [q(end-1,:,:); q(end,:,:); q; q(1,:,:)];
    end
    % qminus(i) = q^-_{i-3/2}^3
    qminus = (qExtended(1:end-1,:)*phiVectorMinus').^3;
    % qplus(i) = q^+_{i-3/2}^3
    qplus  = (qExtended(2:end,:)*phiVectorPlus').^3;
    % qhat(i) = (q^-_{i-1/2}^3 + q^+_{i-1/2}^3)/2
    qhat   = (qminus+qplus)/2;

    E = @(i) 1/deltaX*(getBMatrix(nBasisCpts, q, (i>1)*(i-1) + (i==1)*(strcmp(bc, 'extrapolation') + strcmp(bc, 'periodic')*nCells)) - qminus(i+1)*Phi0 + (qplus(i) - qhat(i))*Phi1);
    F = @(i) 1/deltaX*qhat(i+1)*Phi2;

    U = sparse((nCells+1)*nBasisCpts, (nCells+2)*nBasisCpts);
    for i = 1:nCells+1
        U((nBasisCpts*(i-1)+1):nBasisCpts*i,(nBasisCpts*(i-1)+1):nBasisCpts*i) = E(i);
        U((nBasisCpts*(i-1)+1):nBasisCpts*i,(nBasisCpts*(i)+1):nBasisCpts*(i+1)) = F(i);
    end
end
