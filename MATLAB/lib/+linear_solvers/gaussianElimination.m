function [ x ] = gaussianElimination( A, rhs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [nRows, nCols] = size(A);
    
    for j = 1:nCols
        for i = (j+1):nRows
           scalingConstant = A(i,j)/A(j,j);
           A(i,:) = A(i,:) - scalingConstant*A(j,:);
           rhs(i) = rhs(i) - scalingConstant*rhs(j);
        end
    end

    
        for i = j:-1:1
            for j = nCols:-1:(i+1)
                rhs(i) = rhs(i) - A(i, j)*x(j);
            end
            x(i) = rhs(i)/A(i,i);
        end
            
end

