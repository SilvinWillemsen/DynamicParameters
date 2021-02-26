clear all;
close all;
fs = 44100;
k = 1/fs;

c = 1470;
h = c * k;

N = floor(1/h);
h = 1/N;
N = N-1; % N-1 moving points for N intervals

Dxx = 1/h^2 * (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
sparse(1:N-1, 2:N, ones(1, N-1), N, N));

Q = [2 * eye(N) + c^2 * k^2 * full(Dxx),     -eye(N);
     eye(N),                                zeros(N)];
 
analyseQ (Q, k);

% %% other way
% omega = 2/k * asin(k * c * sqrt(-eig(full(Dxx))));
% f = imag(omega) / (2 * pi);
% scatter(1:length(f), f)

%% Backwards damping
sig0 = 1;
A = eye(N);
B = (2 - 2 * sig0 * k) * eye(N) + c^2 * k^2 * full(Dxx);
C = -(1- 2 * sig0*k) * eye(N);
Q = [A\B,     A\C;
     eye(N),  zeros(N)];

analyseQ (Q, k);

%% Centered damping
sig0 = 1;
A = (1+sig0*k) * eye(N);
B = 2 * eye(N) + c^2 * k^2 * full(Dxx);
C = -(1-sig0*k) * eye(N);
Q = [A\B,     A\C;
     eye(N),  zeros(N)];

analyseQ (Q, k);

%% Adding stiffness
%% Frequency dependent damping (Explicit)
c = 1470;
sig0 = 5;
sig1 = 1;
kappa = 0;
h = sqrt((c^2 * k^2 + 4 * sig1 * k + sqrt((c^2 * k^2 + 4 * sig1 * k)^2 + 16 * kappa^2 * k^2)) / 2);

N = floor(1/h);
h = 1/N;
N = N-1; % N-1 moving points for N intervals

Dxx = 1/h^2 * (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
sparse(1:N-1, 2:N, ones(1, N-1), N, N));
Dxxxx = Dxx * Dxx;

A = (1+sig0*k) * eye(N);
B = 2 * eye(N) + c^2 * k^2 * full(Dxx) - kappa^2 * k^2 * full(Dxxxx) + 2 * sig1 * k * full(Dxx);
C = -(1-sig0*k) * eye(N) - 2 * sig1 * k * Dxx;

Q = [A\B,     A\C;
     eye(N),  zeros(N)];

analyseQ (Q, k)

%% Frequency dependent damping (Implicit)
% sig0 = 5;
% sig1 = 1;
% kappa = 0;
h = sqrt((c^2 * k^2 + sqrt(c^4 * k^4 + 16 * kappa^2 * k^2)) / 2);

N = floor(1/h);
h = 1/N;
N = N-1; % N-1 moving points for N intervals

Dxx = 1/h^2 * (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
sparse(1:N-1, 2:N, ones(1, N-1), N, N));
Dxxxx = Dxx * Dxx;

A = (1+sig0*k) * eye(N) - sig1 * k * full(Dxx);
B = 2 * eye(N) + c^2 * k^2 * full(Dxx) - kappa^2 * k^2 * full(Dxxxx);
C = -(1-sig0*k) * eye(N) - sig1 * k * full(Dxx);

Q = [A\B,     A\C;
     eye(N),  zeros(N)];
analyseQ (Q, k, true)




