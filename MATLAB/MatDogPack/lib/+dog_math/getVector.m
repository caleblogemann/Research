function [vector] = getVector(q)
    [nCells, nBasisCpts] = size(q);
    vector = zeros(nCells*nBasisCpts, 1);
    for i = 1:nCells
        for k = 1:nBasisCpts
            vector(nBasisCpts*(i-1) + k) = q(i, k);
        end
    end
end
