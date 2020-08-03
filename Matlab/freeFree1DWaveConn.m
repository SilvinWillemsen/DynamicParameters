clear all;
close all;
clc;

drawSpeed = 10000;
fs = 44100;
k = 1/fs;
lengthSound = fs;

N = 29.75;
L = 1;
h = 1/N;

cInit = h/k;
c = cInit;

h = c*k;
N = floor(1/h);
h = 1/N;
lambdaSq = (c*k/h)^2

uLeftNext = zeros(floor(N/2) + 1, 1);
uLeft = zeros(floor(N/2) + 1, 1);

if mod(N,2) == 1
    uLeftNext = [uLeftNext; 0];
    uLeft = [uLeft; 0];
end

uLeft(floor(N/4)-4:floor(N/4)+4) = hann(9);
uLeftPrev = uLeft;

uRightNext = zeros(floor(N/2) + 1, 1);
uRight = zeros(floor(N/2) + 1, 1);
% u3(floor(N/4)-5:floor(N/4)+5) = hann(11);
uRightPrev = uRight;

halfRange = 2:length(uLeft)-1;

origHLocs = 0:h:1;

flag = false;
changeC = false;
for n = 1:lengthSound
    NPrev = N;
    if changeC
        c = cInit * (1-0.25*n/lengthSound);
    else
        c = c;
    end
    h = c*k;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    hLocs = 1:-h:0;
    
    lambdaSq = c^2 * k^2 / h^2;

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        if mod(N,2) == 1
            uLeftNext = [uLeftNext; uRightNext(1)];
            uLeft = [uLeft; uRight(1)];
            uLeftPrev = [uLeftPrev; uRightPrev(1)];
        else 
            uRightNext = [uLeftNext(end);uRightNext];
            uRight = [uLeft(end);uRight];
            uRightPrev = [uLeftPrev(end);uRightPrev];
        end
    end
    
    if N < NPrev
        if mod(N,2) == 0
            uLeftNext = uLeftNext(1:end-1);
            uLeft = uLeft(1:end-1);
            uLeftPrev = uLeftPrev(1:end-1);
        else 
            uRightNext = uRightNext(2:end);
            uRight = uRight(2:end);
            uRightPrev = uRightPrev(2:end);
        end
    end

    
    leftRange = 2:length(uLeft)-1;
    rightRange = 2:length(uRight)-1;
    if mod(N,2) == 1
        hLocsLeft = 0:h:((L + h) /2);
        hLocsRight = flip(L:-h:(L/2));
    else
        hLocsLeft = 0:h:(L/2);
        hLocsRight = flip(L:-h:((L - h * 0.9999999)/2));
    end
    alf = (hLocsRight(1) - hLocsLeft(end)) / h;
    uLeftPoint = alf * uLeft(end) + (1 - alf) * uLeft(end-1);
    uRightPoint = alf * uRight(1) + (1 - alf) * uRight(2);

    %% left half string
    uLeftNext(leftRange) = (2-2*lambdaSq) * uLeft(leftRange) + lambdaSq * (uLeft(leftRange+1) + uLeft(leftRange-1)) - uLeftPrev(leftRange);
%     uLeftNext(end) = (2-2*lambdaSq) * uLeft(end) + lambdaSq * (2 * uLeft(end-1)) - uLeftPrev(end);
    uLeftNext(end) = (2-2*lambdaSq) * uLeft(end) + lambdaSq * (uLeft(end-1) + uRightPoint) - uLeftPrev(end);
    
    %% right half string
    uRightNext(rightRange) = (2-2*lambdaSq) * uRight(rightRange) + lambdaSq * (uRight(rightRange+1) + uRight(rightRange-1)) - uRightPrev(rightRange);
%     uRightNext(1) = (2-2*lambdaSq) * uRight(1) + lambdaSq * (2 * uRight(2)) - uRightPrev(1);
    uRightNext(1) = (2-2*lambdaSq) * uRight(1) + lambdaSq * (uRight(2) + uLeftPoint) - uRightPrev(1);

    %% connection
%     F = h * lambdaSq * (uRight(2) - uLeft(end-1));
%     
%     uLeftNext(end) = uLeftNext(end) + 1/h * F;
%     uRightNext(1) = uRightNext(1) - 1/h * F;

    if mod(n, drawSpeed) == 0

        
        hold off;
        plot(hLocsLeft * N, uLeft, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 
%         plot([hLocsLeft, hLocsLeft(end) + h] * N, [0;0;0;0;0;0;-5.55111512312578e-17;-0.146446609406726;-0.353553390593274;-0.500000000000000;-0.500000000000000;-0.353553390593274;-0.146446609406726;0;0;0;0], 'LineWidth' , 2, 'Marker', '.', 'MarkerSize', 20)
        hold on;
        plot(hLocsRight * N, uRight, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
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
        if n == 16
            disp("wait")
        end
        drawnow;
    end
    uLeftPrev = uLeft;
    uLeft = uLeftNext;
    
    uRightPrev = uRight;
    uRight = uRightNext;
    
    outFree(n) = uRightNext(end - 10);
end
subplot(2,1,1)
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])

c/2