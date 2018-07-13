function [qVector] = getQVector(q)
    [nCells, nBasisCpts] = size(q);
    qVector = zeros(nCells*nBasisCpts, 1);
    for i = 1:nCells
        for k = 1:nBasisCpts
            qVector(nBasisCpts*(i-1) + k) = q(i, k);
        end
    end
end
