clear all;
close all;
clc;

drawSpeed = 1;
fs = 44100;
k = 1/fs;
lengthSound = fs;

N = 30;
L = 1;
h = 1/N;

c = h/k;

h = c*k;
N = floor(1/h);
h = 1/N;
lambdaSq = (c*k/h)^2


uNext = zeros((N - 1) * 2, 1);
u = zeros((N - 1) * 2, 1);

u(floor(N/4)-4:floor(N/4)+4) = hann(9);
uPrev = u;

e = ones(length(u), 1);
DxxPre = spdiags([e -2*e e], -1:1, N-1, N-1);

Dxx = zeros(length(u), length(u));
Dxx(1:N-1, 1:N-1) = DxxPre;
Dxx(N:end, N:end) = DxxPre;

for n = 1:lengthSound

    alpha = 0.1;
    cubicInterp = [alpha * (alpha-1) * (alpha-2) / -6; ...
        (alpha - 1) * (alpha + 1) * (alpha - 2) / 2; ...
        alpha * (alpha + 1) * (alpha - 2) / -2; ...
        alpha * (alpha + 1) * (alpha - 1) / 6 ];
    I = zeros(1, N-1);  
    I(floor((N-1) / 2) - 2 : floor((N-1) / 2) + 1) = cubicInterp';

    J = zeros((N-1) * 2, 1);
    J(1:N-1) = I' * 1/h;
    J(N:end) = -I' * 1/h;

    uNext = 2 * u + lambdaSq * Dxx * u - uPrev;
   
    uNext = uNext + J .* F;
    %% connection
%     F1 = h/2 * lambdaSq * (Iw1 * Dxxw * w - I * Dxxu * u);
%     F2 = h/2 * lambdaSq * (Iw2 * Dxxw * w - Iu2 * Dxxu * u);
%     F2 = 0;
%     F = h * lambdaSq * (uRight(2) - uLeft(end-1));
 
%     uNext = uNext + (J * F1);% + Ju2 * F2);
%     wNext = wNext - (Jw1 * F1);% + Jw2 * F2);

    if mod(n, drawSpeed) == 0

        
        hold off;
        plot(u(1:N-1), 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 
%         plot([hLocsLeft, hLocsLeft(end) + h] * N, [0;0;0;0;0;0;-5.55111512312578e-17;-0.146446609406726;-0.353553390593274;-0.500000000000000;-0.500000000000000;-0.353553390593274;-0.146446609406726;0;0;0;0], 'LineWidth' , 2, 'Marker', '.', 'MarkerSize', 20)
        hold on;
        plot(u(N:end), 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
%         plot(hLocsRight * N, [0;5.55111512312578e-17;0;0;0.146446609406726;0.353553390593274;0.500000000000000;0.500000000000000;0.353553390593274;0.146446609406726;0;0;0;0;0;0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)
        ylim([-0.6, 0.6])
        grid on;
%         hold on;
%         scatter([hLocsLeft(end) + h, hLocsRight(1) - h], [uRightPoint, uLeftPoint]);
%         scatter((hLocsLeft(end)+h) * N, 0, 'b', 'Marker', 'o', 'Linewidth', 2)
%         scatter((hLocsRight(1)-h) * N, 0, 'r', 'filled')
%         legend('$u$', '$w$', '$u_{M+1}$', '$w_{-1}$', 'interpreter', 'latex', 'Fontsize', 16)
        legend('$u$', '$w$', 'interpreter', 'latex', 'Fontsize', 16)

        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 6
            disp("wait")
        end
        drawnow;
    end
    uPrev = u;
    u = uNext;
    
    outFree(n) = uNext(end - 10);
end
subplot(2,1,1)
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])

c/2