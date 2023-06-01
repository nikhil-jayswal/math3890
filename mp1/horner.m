function v = horner(c, t)
    % Nikhil Jayswal
    % MATH 3890
    % Machine Problem 1, Question 1
    % 2/1/2021
    
    % c = set of (n+1) coeffs.
    % t = set of points where polynomial is to be evaluated 
    
    v = zeros(size(t)); % # of output values = # of evaluation points
    n = length(c);
    for i = 1:length(t)
        x = t(i);
        for j = 1:length(c)
            v(i) = x*v(i) + c(n - j + 1);
        end
    end
end