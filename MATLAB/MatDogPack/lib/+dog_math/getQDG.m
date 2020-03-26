function [qDG] = getQDG(qFD)
    qDG = zeros(length(qFD), 1, 2);
    qDG(:, 1, 1) = qFD(:);
    qDG(:, 1, 2) = (qFD([2:end, 1]) - qFD([end,1:end-1]))*(sqrt(3)/12);
end
