function [c, M] = scatrbf(x, y, z, eps, rbf)
    
    % number of basis functions/coefficients = length of x/y/z
    n =length(z);
    
    % set up observation equations
    M = ones(n, n);
    for i = 1:n
        for j = 1:n
            % M(i, j) = 1 when i = j; don't compute
            if i ~= j
                % compute distance
                r = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
                M(i, j) = rbf(eps, r);
            end
        end
    end
    
    % compute coefficients
    c = M\z;
end