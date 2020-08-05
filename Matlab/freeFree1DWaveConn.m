clear all;
close all;
clc;

drawSpeed = 1;
fs = 44100;
k = 1/fs;
lengthSound = fs;

Ninit = 30;
N = Ninit;
L = 1;
h = 1/N;

cInit = h/k;
c = cInit;

h = c*k;
N = floor(1/h);
h = 1/N;
lambdaSq = (c*k/h)^2

uNext = zeros(N-1, 1);
u = zeros(N-1, 1);

if mod(N,2) == 1
    uNext = [uNext; 0];
    u = [u; 0];
end

u(floor(N/5)-4:floor(N/5)+4) = hann(9);
uPrev = u;

wNext = zeros(N - 1, 1);
w = zeros(N - 1, 1);
% u3(floor(N/4)-5:floor(N/4)+5) = hann(11);
% w(floor(N/5)-4:floor(N/5)+4) = hann(9);

wPrev = w;

origHLocs = 0:h:1;

eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
% Dxxu(end, end - 1) = 2;
ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
% Dxxw(end, end - 1) = 2;

flag = false;
changeC = false;

eta = 0;
etaNext = 0;
omega0 = 1000000;
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
            uNext = [uNext; wNext(1)];
            u = [u; w(1)];
            uPrev = [uPrev; wPrev(1)];
        else 
            wNext = [uNext(end);wNext];
            w = [u(end);w];
            wPrev = [uPrev(end);wPrev];
        end
    end
    
    if N < NPrev
        if mod(N,2) == 0
            uNext = uNext(1:end-1);
            u = u(1:end-1);
            uPrev = uPrev(1:end-1);
        else 
            wNext = wNext(2:end);
            w = w(2:end);
            wPrev = wPrev(2:end);
        end
    end

    
    leftRange = 2:length(u)-1;
    rightRange = 2:length(w)-1;
%     if mod(N,2) == 1
% %         hLocsLeft = 0:h:((L + h * 0.99999) /2);
% %         hLocsRight = flip(L:-h:(L/2));
%         hLocsLeft = 0:h:(L/2 + h * (overlap + 0.9999));
%         hLocsRight = flip(L:-h:(L/2 - h * overlap));
% 
%     else
%         hLocsLeft = 0:h:(L/2 + h * overlap);
%         hLocsRight = flip(L:-h:(L/2 - h * (overlap)));
%     end
    hLocsLeft = (1:(length(u))) * h;
    hLocsRight = flip(1 - ((1:(length(w))) * h));

    alf = 0.25;
    I = zeros(1, length(u));
    I(floor(N/2)+1) = alf;
    I(floor(N/2)) = 1-alf;
    J = I' * 1/h;
        
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;

    %% connection
    etaPrev = eta;
    eta = I * u - I * w;
%     eta - etaNext
    F = ((I * wNext - I * uNext) - etaPrev) / ((I * J + I * J) * k^2 + 2/omega0^2);        

%     etaNext = (2 * eta - etaPrev - k^2 * omega0^2 / h * etaPrev + (I * uNext - I * wNext)) / (1 + k^2 * omega0^2/h);
%     F = (-eta - etaPrev) / (2 * k^2 + 2/omega0^2);
%     F = -omega0^2 * (etaNext + etaPrev) / 2;
%     Ftest = h * ((c^2  / h^2 * Iw * Dxxw * w - c^2  / h^2 * I * Dxxu * u) + 1/k^2 * (etaNext - 2 * eta + etaPrev)) / 2;
%     F - Ftest
%     F = h / 2 * c^2 / (h^2) * (Iw * Dxxw * w - I * Dxxu * u);
    uNext = uNext + k^2 * J * F;
    wNext = wNext - k^2 * J * F;
    trueEtaSave(n) = (I * uNext - I * wNext);
    etaNext = -2 / omega0^2 - etaPrev;
    etaSave(n) = etaNext;
    trueEtaSave(n) - etaSave(n)
%     solut = [2/h, 2 * Iu1 * Ju2; 2 * Iw2 * Jw1, 2/h] \ [lambdaSq * (Iw1 * Dxxw * w - Iu1 * Dxxu * u); lambdaSq * (Iw2 * Dxxw * w - Iu2 * Dxxu * u)] ;
%     F1 = h/2 * lambdaSq * (sum(Iw1 * Dxxw * w) - sum(Iu1 * Dxxu * u));
%     F2 = h/2 * lambdaSq * (sum(Iw2 * Dxxw * w) - sum(Iu2 * Dxxu * u));
%     uNext = uNext + (Ju1 * solut(1) + Ju2 * solut(2));
%     wNext = wNext - (Jw1 * solut(1) + Jw2 * solut(2));
    scalingU = ones(length(u),1);
%     scalingU(end) = 0.5;
    kinEnergyU(n) = 1/2 * h * sum(scalingU .* (1/k * (u - uPrev)).^2);
    potEnergyU(n) = c^2/(2 * h) * sum((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)));
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0)); % left boundary
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((0 - u(end)) .* (0 - uPrev(end))); % right boundary

%     potEnergyBoundU(n) = c^2/(2 * h) * sum((u(end-1) - u(end)) .* (uPrev(end-1) - uPrev(end))); % right boundary

    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);% + potEnergyBoundU(n);
    
    scalingW = ones(length(w),1);
%     scalingW(end) = 0.5;

    kinEnergyW(n) = 1/2 * h * sum(scalingW .* (1/k * (w - wPrev)).^2);
    potEnergyW(n) = c^2/(2 * h) * sum((w(2:end) - w(1:end-1)) .* (wPrev(2:end) - wPrev(1:end-1)));
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((w(1) - 0) .* (wPrev(1) - 0)); % left boundary
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % right boundary

%     potEnergyBoundW(n) = c^2/(2 * h) * sum((w(end-1) - w(end)) .* (wPrev(end-1) - wPrev(end))); % right boundary

%     potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % left boundary
%     potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((w(1) - w(2)) .* (wPrev(1) - wPrev(2))); % right boundary

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);% + potEnergyBoundW(n);
    
    connEnergy(n) = omega0^2 / 2 * 1/2 * (eta^2 + etaPrev^2);

    totEnergy(n) = totEnergyU(n) + totEnergyW(n) + connEnergy(n);
    if mod(n, drawSpeed) == 0

        subplot(3,1,1)
        hold off;
        plot(u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 
%         plot([hLocsLeft, hLocsLeft(end) + h] * N, [0;0;0;0;0;0;-5.55111512312578e-17;-0.146446609406726;-0.353553390593274;-0.500000000000000;-0.500000000000000;-0.353553390593274;-0.146446609406726;0;0;0;0], 'LineWidth' , 2, 'Marker', '.', 'MarkerSize', 20)
        hold on;
        plot(w, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
%         plot(hLocsRight * N, [0;5.55111512312578e-17;0;0;0.146446609406726;0.353553390593274;0.500000000000000;0.500000000000000;0.353553390593274;0.146446609406726;0;0;0;0;0;0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)
        ylim([-1, 1])
        grid on;
%         hold on;
%         scatter([hLocsLeft(end) + h, hLocsRight(1) - h], [uRightPoint, uLeftPoint]);
%         scatter((hLocsLeft(end)+h) * N, 0, 'b', 'Marker', 'o', 'Linewidth', 2)
%         scatter((hLocsRight(1)-h) * N, 0, 'r', 'filled')
%         legend('$u$', '$w$', '$u_{M+1}$', '$w_{-1}$', 'interpreter', 'latex', 'Fontsize', 16)
%         legend('$u$', '$w$', 'interpreter', 'latex', 'Fontsize', 16)

        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 6
            disp("wait")
        end
        
        subplot(3,1,2)
        hold off
%         plot(potEnergyU(1:n) + kinEnergyU(1:n) + potEnergyW(1:n) + kinEnergyW(1:n) - (potEnergyU(1) + kinEnergyU(1) + potEnergyW(1) + kinEnergyW(1)));
%         hold on;
        plot(totEnergyU(1:n) + totEnergyW(1:n) - (totEnergyU(1) + totEnergyW(1)))   
%         plot(potEnergyBoundW(1:n))     

        hold on;
        plot(connEnergy(1:n))
%         hold off;
        
%         hold off;
%         plot(totEnergyU(1:n) / totEnergy(1) - 1);
%         hold on;
%         plot(connEnergy(1:n) / 2)
        subplot(3,1,3);
        plot(totEnergy(1:n) - totEnergy(1));

        drawnow;
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
    outFree(n) = wNext(end - 10);
end
subplot(2,1,1)
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])

c/2