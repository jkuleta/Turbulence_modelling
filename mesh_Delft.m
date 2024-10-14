function [MESH] = mesh_Delft(H, n)

% Parameters
fact = 6.0; 
ns = 3; 

% Finite Difference Method
di = 1.0 / (n - 1);
i = (0:n-1)' / (n - 1);  % Normalized coordinates from 0 to 1

% y - coordinate: tanh clustering at the walls
y = H * (1.0 + tanh(fact * (i - 0.5)) / tanh(fact / 2)) / 2.0;
y = (y - y(1)) / (y(end) - y(1));  % Normalize y to range from 0 to 1

% Coordinate transformation: derivative of y with respect to 'i'
dydi = H * fact / 2.0 / tanh(fact / 2) ./ cosh(fact * (i - 0.5)).^2.0;
dydi = (dydi + dydi(end:-1:1)) / 2.0;  % Adjust for symmetry

% Coordinate transformation: second derivative of y with respect to 'i'
d2ydi2 = -H * fact^2.0 / tanh(fact / 2) .* tanh(fact * (i - 0.5)) ./ cosh(fact * (i - 0.5)).^2.0;
d2ydi2 = (d2ydi2 - d2ydi2(end:-1:1)) / 2.0;  % Adjust for symmetry

% Coefficient matrix for d()/dy
ddy = zeros(n, n);
ddy(1, 1:7) = finiteDiffCoeff(0:6, 1);
ddy(2, 1:7) = finiteDiffCoeff(-1:5, 1);
ddy(3, 1:7) = finiteDiffCoeff(-2:4, 1);
ddy(n-2, n-6:n) = finiteDiffCoeff(-4:2, 1);
ddy(n-1, n-6:n) = finiteDiffCoeff(-5:1, 1);
ddy(n, n-6:n) = finiteDiffCoeff(-6:0, 1);

% Central finite difference for interior points
for i = 1+ns:n-ns
    ddy(i, :) = 0.0;
    ddy(i, i-ns:i+ns) = finiteDiffCoeff(-ns:ns, 1);
end

% Multiply coordinate transformation
ddy = bsxfun(@times, ddy, 1 / di ./ dydi);

% Coefficient matrix for d2()/dy2 (second derivative)
d2dy2 = zeros(n, n);
d2dy2(1, 1:7) = finiteDiffCoeff(0:6, 2);
d2dy2(2, 1:7) = finiteDiffCoeff(-1:5, 2);
d2dy2(3, 1:7) = finiteDiffCoeff(-2:4, 2);
d2dy2(n-2, n-6:n) = finiteDiffCoeff(-4:2, 2);
d2dy2(n-1, n-6:n) = finiteDiffCoeff(-5:1, 2);
d2dy2(n, n-6:n) = finiteDiffCoeff(-6:0, 2);

% Central finite difference for second derivatives
for i = 1+ns:n-ns
    d2dy2(i, :) = 0.0;
    d2dy2(i, i-ns:i+ns) = finiteDiffCoeff(-ns:ns, 2);
end

% Multiply coordinate transformation for second derivatives
d2dy2 = bsxfun(@times, d2dy2, 1 / di^2 ./ dydi.^2) ...
        - bsxfun(@times, ddy, d2ydi2 ./ dydi.^2);

% Create mesh structure
MESH = struct('y', y, 'ddy', ddy, 'd2dy2', d2dy2);

end
