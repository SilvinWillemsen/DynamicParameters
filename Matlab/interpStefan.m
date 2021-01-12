% optimised interpolant over unequally spaced grid points

% clear all
hold off

%close all

N = 11;             % number of points
alpha = 1;        % relative bandwidth range
freq = 5;          % spatial frequency

x = rand(N,1)-0.5;                  % randomly select grid points
% x = -0.52 + [0:N-2]'/(N-2);           % choose this regularly spacing...
% x = [x(1:floor(N/2)); x(floor(N/2)); x(floor(N/2) + 1 : end)];
% x = xRand;
x = sort(x);
h = max(diff(x));

y = sin(freq*x);                    % use a sine wave for the function

bmax = alpha*pi/h;

b = (sin(bmax*x)./x);
dist = x*ones(1,N)-ones(N,1)*x';    % distance matrix between points
A = sin(bmax*dist)./dist;
A(1+(N+1)*[0:N-1]') = bmax;         % collection of sinc functions with centers at grid point locations

a = A\b;           % optimal coefficients
y_est = a'*y       % form interpolant

xplot = [-0.5:0.01:0.5];
plot(x,y,'kx', xplot, sin(freq*xplot), 'g')
hold on
plot(0,y_est,'rx')
y_est
