function [phi] = legendrePolynomialDerivative(k, xi)
    if(k == 1)
        phi = 0;
    elseif (k==2)
        phi = sqrt(3);
    elseif (k==3)
        phi = sqrt(5)*3*xi;
    elseif (k==4)
        phi = sqrt(7)/2*(15*xi.^2 - 3);
    else
        error('This order polynomial is not implemented');
    end
end
