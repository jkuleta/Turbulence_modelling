function [MESH] = mesh(H, n)
    % Generate uniform mesh (n points between 0 and H)
    y = linspace(0, H, n)';  % n equally spaced points from 0 to H
    di = y(2) - y(1);  % uniform spacing between nodes

    % First derivative matrix ddy using central differences
    ddy = zeros(n, n);
    ddy(1, 1:3) = finiteDiffCoeff(0:2, 1);  % Forward difference at the first point
    ddy(n, n-2:n) = finiteDiffCoeff(-2:0, 1);  % Backward difference at the last point
          
    for i = 2:n-1
    ddy(i, i-1:i+1) = finiteDiffCoeff(-1:1, 1);  % Central difference for internal points
    end
    ddy = ddy / di;  

    % Second derivative matrix d2dy2 using central differences
    d2dy2 = zeros(n, n);
    d2dy2(1, 1:3) = finiteDiffCoeff(0:2, 2);  % Forward difference at the first point
    d2dy2(n, n-2:n) = finiteDiffCoeff(-2:0, 2);  % Backward difference at the last point
           
    for i = 2:n-1
    d2dy2(i, i-1:i+1) = finiteDiffCoeff(-1:1, 2);  % Central difference for internal points
    end
    d2dy2 = d2dy2 / di^2; 
    
    MESH = struct('y', y, 'ddy', ddy, 'd2dy2', d2dy2);

end

