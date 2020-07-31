clear all;
close all;

fs = 44100;
k = 1/fs;
lengthSound = fs;

c = 1000;
h = c*k;

L = 1;
N = floor (L/h);
h = 1/N;

N = N-1;

lambdaSq = (c * k / h)^2

uNext = zeros (N,1);
u = zeros (N,1);
u(floor(N/2 - 3): floor(N/2 + 3)) = hann(7);
uPrev = u;

e = ones(N,1);
Dxx = spdiags([e -2*e e], -1:1,N,N);
alpha = 0.5;
for n = 1:lengthSound
    uNext = 2 * u + lambdaSq * Dxx * u - uPrev;
    
    x = [-h, 0, h*(1+alpha)];    
    hold off;
    plot([(0:N) * h, 1], [0; uNext; 0], 'Marker', '.', 'MarkerSize', 20);
    hold on;
    plot([(N - 1)*h,N*h, 1], uNext(end) + 1/(h * (alpha + 2)) * (0 - uNext(end-1)) * x + (0 - (alpha + 2) * uNext(end) + (alpha + 1) * uNext(end-1))/(h^2 * (alpha + 2) * (alpha + 1)) * x.^2)  
    xlim([1-4*h, 1]);
    drawnow;
    pause(0.2)

    uPrev = u;
    u = uNext;
end
