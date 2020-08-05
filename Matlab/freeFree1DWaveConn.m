clear all;
close all;
clc;

drawSpeed = 10000;
fs = 44100;
k = 1/fs;
lengthSound = fs;

Ninit = 31.5;
N = Ninit;
L = 1;
h = 1/N;

cInit = h/k;
c = cInit;

h = c*k;
N = floor(1/h);
h = 1/N;
lambdaSq = (c*k/h)^2

overlap = 2;

uNext = zeros(floor(N/2) + ceil((overlap + 1) / 2) - 1, 1);
u = zeros(floor(N/2) + ceil((overlap + 1) / 2) - 1, 1);

if mod(N,2) == 1
    uNext = [uNext; 0];
    u = [u; 0];
end

u(floor(N/5)-4:floor(N/5)+4) = hann(9);
uPrev = u;

wNext = zeros(floor(N/2) + floor((overlap + 1) / 2) - 1, 1);
w = zeros(floor(N/2) + floor((overlap + 1) / 2) - 1, 1);
% u3(floor(N/4)-5:floor(N/4)+5) = hann(11);
% w(floor(2*N/5)-4:floor(2*N/5)+4) = hann(9);

wPrev = w;

origHLocs = 0:h:1;

eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
Dxxu(end, end - 1) = 2;
ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
Dxxw(1, 2) = 2;

flag = false;
changeC = false;

eta = 0;
etaNext = 0;
omega0 = 10000;
for n = 1:lengthSound
    NPrev = N;
    if changeC
        c = cInit * (1-5*n/lengthSound);
    else
        c = c;
    end
    h = c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    hLocs = 1:-h:0;
    
    lambdaSq = c^2 * k^2 / h^2;
    hLocsLeft = (1:(length(u))) * h;
    hLocsRight = flip(1 - ((1:(length(w))) * h));

%     alf = (Ninit - N) + 0.5 * (1-(Ninit - N));
    alf = 1 - 0.5 * (hLocsLeft(end) - hLocsRight(1)) / h;
    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        if mod(N,2) == 1
            uNext = [uNext; ((1-alf) * wNext(1) + alf * wNext(2) - alf * uNext(end)) / (1-alf)];
            u = [u; ((1-alf) * w(1) + alf * w(2) - alf * u(end)) / (1-alf)];
            uPrev = [uPrev; ((1-alf) * wPrev(1) + alf * wPrev(2) - alf * uPrev(end)) / (1-alf)];
        else 
            wNext = [uNext(end - overlap + 1);wNext];
            w = [u(end - overlap + 1);w];
            wPrev = [uPrev(end - overlap + 1);wPrev];
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        Dxxu(end, end - 1) = 2;
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
        Dxxw(1, 2) = 2;
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
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        Dxxu(end, end - 1) = 2;
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
        Dxxw(1, 2) = 2;
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

%     alf = (hLocsRight(1 + overlap) - hLocsLeft(end - overlap)) / h;
%     alf = 0.25;
    Iu = zeros(1, length(u));
    Iu(end) = alf;
    Iu(end-1) = 1-alf;
    Ju = Iu' * 1/h;
%     alf = 1;
    Iw = zeros(1, length(w));
    Iw(2) = 1-alf;
    Iw(1) = alf;
    Jw = Iw' * 1/h;

        
    uLeftLocs = ((length(uNext) - overlap + 1) : length(uNext))';
    uRightLocs = (1:overlap)';
    
    uLeftPoints = alf * u(uLeftLocs) + (1 - alf) * u(uLeftLocs - 1);
    uRightPoints = alf * w(uRightLocs) + (1 - alf) * w(uRightLocs + 1);
    
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;

    %% connection
%     etaPrev = eta;
%     eta = Iu * u - Iw * w;
% %     eta - etaNext
%     etaNext = (2 * eta - etaPrev - k^2 * omega0^2 / h * etaPrev + k^2 * (c^2 / h^2 * Iu * Dxxu * u - c^2 / h^2 * Iw * Dxxw * w)) / (1 + k^2 * omega0^2/h);
%     F = -omega0^2 * (etaNext + etaPrev) / 2;
%     Ftest = h * ((c^2  / h^2 * Iw * Dxxw * w - c^2  / h^2 * Iu * Dxxu * u) + 1/k^2 * (etaNext - 2 * eta + etaPrev)) / 2;
%     F - Ftest
    F = c^2 / (h^2) * (Iw * Dxxw * w - Iu * Dxxu * u) / (Iu * Ju + Iw * Jw);
    uNext = uNext + k^2 * Ju * F;
    wNext = wNext - k^2 * Jw * F;
    
%     eta = Iu * u - Iw * w
    
%     trueEtaSave(n) = (Iu * uNext - Iw * wNext);
%     etaSave(n) = etaNext;
%     solut = [2/h, 2 * Iu1 * Ju2; 2 * Iw2 * Jw1, 2/h] \ [lambdaSq * (Iw1 * Dxxw * w - Iu1 * Dxxu * u); lambdaSq * (Iw2 * Dxxw * w - Iu2 * Dxxu * u)] ;
%     F1 = h/2 * lambdaSq * (sum(Iw1 * Dxxw * w) - sum(Iu1 * Dxxu * u));
%     F2 = h/2 * lambdaSq * (sum(Iw2 * Dxxw * w) - sum(Iu2 * Dxxu * u));
%     uNext = uNext + (Ju1 * solut(1) + Ju2 * solut(2));
%     wNext = wNext - (Jw1 * solut(1) + Jw2 * solut(2));
    scalingU = ones(length(u),1);
    scalingU(end) = 0.5;
    kinEnergyU(n) = 1/2 * h * sum(scalingU .* (1/k * (u - uPrev)).^2);
    potEnergyU(n) = c^2/(2 * h) * sum((u(2:end-1) - u(1:end-2)) .* (uPrev(2:end-1) - uPrev(1:end-2)));
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0)); % left boundary
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(end-1) - u(end)) .* (uPrev(end-1) - uPrev(end))); % right boundary

    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    scalingW = ones(length(w),1);
    scalingW(1) = 0.5;

    kinEnergyW(n) = 1/2 * h * sum(scalingW .* (1/k * (w - wPrev)).^2);
    potEnergyW(n) = c^2/(2 * h) * sum((w(3:end) - w(2:end-1)) .* (wPrev(3:end) - wPrev(2:end-1)));
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % left boundary
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((w(1) - w(2)) .* (wPrev(1) - wPrev(2))); % right boundary

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);
    
%     connEnergy(n) = omega0^2 / 2 * 1/2 * (eta^2 + etaPrev^2);
    
    totEnergy(n) = totEnergyU(n) + totEnergyW(n);% + connEnergy(n);
    if mod(n, drawSpeed) == 0

%         subplot(3,1,1)
        hold off;
        plot(hLocsLeft * N, u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 
%         plot([hLocsLeft, hLocsLeft(end) + h] * N, [0;0;0;0;0;0;-5.55111512312578e-17;-0.146446609406726;-0.353553390593274;-0.500000000000000;-0.500000000000000;-0.353553390593274;-0.146446609406726;0;0;0;0], 'LineWidth' , 2, 'Marker', '.', 'MarkerSize', 20)
        hold on;
        plot(hLocsRight * N, w, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
%         plot(hLocsRight * N, [0;5.55111512312578e-17;0;0;0.146446609406726;0.353553390593274;0.500000000000000;0.500000000000000;0.353553390593274;0.146446609406726;0;0;0;0;0;0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)
        if n == 10
            disp("wait");
        end
        ylim([-1, 1])
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
        
%         subplot(3,1,2)
%         hold off
%         plot(totEnergyU(1:n) + totEnergyW(1:n) - (totEnergyU(1) + totEnergyW(1)))     
% %         hold on
% %         plot(connEnergy(1:n))
% %         hold off;
%         
% %         hold off;
% %         plot(totEnergyU(1:n) / totEnergy(1) - 1);
% %         hold on;
% %         plot(connEnergy(1:n) / 2)
%         subplot(3,1,3);
%         plot(totEnergy(1:n) / totEnergy(1) - 1);

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