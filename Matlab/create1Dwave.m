function [B, C, N, h, Dxx] = create1Dwave(c, L, k)
   
    h = c*k;
    N = floor(L/h); % Number of gridpoints
    h = L/N; % Recalculate gridspacing
    
    N = N - 2;

    Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
    Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%     Dxxxx = Dxx*Dxx;
    N = N - 2;
    B = 2 * eye(N) + c^2 * k^2 / h^2 * Dxx(2:end-1,2:end-1);
    C = -1 * eye(N);
    
    N = N + 4;
    
    Dxx2 = zeros(N);
    Dxx2(2:end-1, 2:end-1) = Dxx;
    Dxx = Dxx2;
    
    Dxxxx2 = zeros(N);
    Dxxxx2(2:end-1, 2:end-1) = Dxxxx;
    Dxxxx = Dxxxx2;
    
    B2 = zeros(N);
    B2(3:end-2, 3:end-2) = B;
    B = B2;

    C2 = zeros(N);
    C2(3:end-2, 3:end-2) = C;
    C = C2;
end