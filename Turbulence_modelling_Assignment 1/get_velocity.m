function [u] = get_velocity(DP_dx, visc_nu, visc_nu_eddy, rho, dy, y)
    N = length(y);
    U = ones(N,1);
    A = zeros(N,N);
    B = ones(N,1);
    
    B = ((1./((visc_nu+visc_nu_eddy)*rho)).*DP_dx.*dy^2) .*B;
    %Boundary conditions
    B(1) = 0;
    B(end) = 0;
    
    for i = 2:N-1
        A(i, i-1) = 1;    % Coefficient U(y-dy)
        A(i, i) = -2;       % Coefficient U(y)
        A(i, i+1) = 1;     % Coefficient U(y+dy)
    end
    
    A(1, 1) = 1;
    A(end, end-1) = -1;
    A(end, end) = 1;
    
    u = linsolve(A,B);

end
