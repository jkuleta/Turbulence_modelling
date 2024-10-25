function [MESH] = meshGenerator(H, n)

ns = 3; 

    % y - coordinate: uniformly spaced points between 0 and H
    y = linspace(0, H, n)';
    
    % constant spacing
    dy = y(2) - y(1);
    
    % -------------------------------------------------------------
    % coefficient matrix for d()/dy (first derivative)
    ddy = zeros(n, n);
    
    % Finite difference for first derivative
    ddy(1,1:7) = finiteDiffCoeff(0:6, 1);        % near the boundaries
    ddy(2,1:7) = finiteDiffCoeff(-1:5, 1);
    ddy(3,1:7) = finiteDiffCoeff(-2:4, 1);
    ddy(n-2,n-6:n) = finiteDiffCoeff(-4:2, 1);
    ddy(n-1,n-6:n) = finiteDiffCoeff(-5:1, 1);
    ddy(n,  n-6:n) = finiteDiffCoeff(-6:0, 1);  
    
    for i=1+ns:n-ns
        ddy(i,:) = 0.0;
        ddy(i,i-ns:i+ns) = finiteDiffCoeff(-ns:ns, 1);
    end
    
    ddy = ddy / dy; % Adjust by constant spacing (now dense matrix)
    
    % -------------------------------------------------------------
    % coefficient matrix for d2()/dy2 (second derivative)
    d2dy2 = zeros(n, n);
    
    % Finite difference for second derivative
    d2dy2(1,1:7) = finiteDiffCoeff(0:6, 2);      % near the boundaries
    d2dy2(2,1:7) = finiteDiffCoeff(-1:5, 2);
    d2dy2(3,1:7) = finiteDiffCoeff(-2:4, 2);
    d2dy2(n-2,n-6:n) = finiteDiffCoeff(-4:2, 2);
    d2dy2(n-1,n-6:n) = finiteDiffCoeff(-5:1, 2);
    d2dy2(n,  n-6:n) = finiteDiffCoeff(-6:0, 2);
    
    for i=1+ns:n-ns
        d2dy2(i,:) = 0.0;
        d2dy2(i,i-ns:i+ns) = finiteDiffCoeff(-ns:ns, 2);
    end
    
    d2dy2 = d2dy2 / dy^2; % Adjust by constant spacing (now dense matrix)
    
    % Return the mesh structure
    MESH = struct('y', y, 'ddy', ddy, 'd2dy2', d2dy2);

end
