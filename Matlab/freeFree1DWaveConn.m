clear all;
close all;
clc;

drawSpeed = 10000;
fs = 44100;
k = 1/fs;
lengthSound = fs;

Ninit = 30.5;
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

zNext = zeros(N, 1);
z = zeros(N-1, 1);
z(floor(N/5)-4:floor(N/5)+4) = hann(9);
zPrev = z;

origHLocs = 0:h:1;

eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

Dxxu(end-1,end-2) = 1.5;
Dxxu(end,end-1) = 1.5;


ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
% Dxxw(3,2) = 0.5;
Dxxw(2,3) = 1.5;
Dxxw(1,2) = 1.5;
% Dxxw(2,1) = 1/3;
% Dxxw(1,1) = -2;

ez = ones(length(z), 1);
Dxxz = spdiags([ez -2*ez ez], -1:1, length(z),length(z));

flag = false;
changeC = false;

eta = 0;
etaNext = 0;
omega0 = 10000;
interpol = "linear";
outFree = zeros(lengthSound, 1);
for n = 1:lengthSound
    NPrev = N;
    if changeC
        c = cInit * (1-5*n/lengthSound);
    else
        c = c;
    end
    h = 1.001*c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    hLocs = 1:-h:0;
    
    lambdaSq = c^2 * k^2 / h^2;
    hLocsLeft = (1:(length(u))) * h;
    hLocsRight = flip(1 - ((1:(length(w))) * h));

    alf = (Ninit - N) + 0.5 * (1-(Ninit - N));
%     alf = 1 - 0.5 * (hLocsLeft(floor(N/2)) - hLocsRight(length(hLocsRight) - floor(N/2) + 1)) / h;
  
    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        if mod(N,2) == 1
            if interpol == "linear"
                uNext = [uNext; ((1-alf) * wNext(1) + alf * wNext(2) - alf * uNext(end)) / (1-alf)];
                u = [u; ((1-alf) * w(1) + alf * w(2) - alf * u(end)) / (1-alf)];
                uPrev = [uPrev; ((1-alf) * wPrev(1) + alf * wPrev(2) - alf * uPrev(end)) / (1-alf)];
            else
                uNext = [uNext; ((1-alf) * wNext(1) + alf * wNext(2) - alf * uNext(end)) / (1-alf)];
                u = [u; ((1-alf) * w(1) + alf * w(2) - alf * u(end)) / (1-alf)];
                uPrev = [uPrev; ((1-alf) * wPrev(1) + alf * wPrev(2) - alf * uPrev(end)) / (1-alf)];
            
            end
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


    cub = [alf * (alf - 1) * (alf - 2) / -6, ...
            (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
            alf * (alf + 1) * (alf - 2) / -2, ...
            alf * (alf + 1) * (alf - 1) / 6];
        
    Iu = zeros(1, length(u));
    Iw = zeros(1, length(w));
    if interpol == "linear"
        Iu(end) = alf;
        Iu(end-1) = 1-alf;
        Iw(2) = 1-alf;
        Iw(1) = alf;
    else
        Iu(end-3 : end) = cub;
        Iw(1:4) = cub;
    end
    Ju = Iu' * 1/h;    
    Jw = Iw' * 1/h;

        
    uLeftLocs = ((length(uNext) - overlap + 1) : length(uNext))';
    uRightLocs = (1:overlap)';
    
    uLeftPoints = alf * u(uLeftLocs) + (1 - alf) * u(uLeftLocs - 1);
    uRightPoints = alf * w(uRightLocs) + (1 - alf) * w(uRightLocs + 1);
    
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    
    %% orig string
    zNext = 2 * z + lambdaSq * Dxxz * z - zPrev;


    alf = 0.5;
    %% connection  
    I1u = zeros(1,length(u));
    I1w = zeros(1,length(w));
    I1u(end-1) = 1-alf;
    I1u(end) = alf;
    I1w(1) = 1;
    J1u = I1u' * 1/h;    
    J1w = I1w' * 1/h;

    I2u = zeros(1,length(u));
    I2w = zeros(1,length(w));
    I2u(end) = 1;
    I2w(1) = alf;
    I2w(2) = 1-alf;
    J2u = I2u' * 1/h;    
    J2w = I2w' * 1/h;

    
    a11 = I1u * J1u + I1w * J1w;
    a12 = I1u * J2u + I1w * J2w;
    a21 = I2u * J1u + I2w * J1w;
    a22 = I2u * J2u + I2w * J2w;
    
%     F1 = h * c^2 / h^2 * (-u(end-2));
%     F2 = h * c^2 / h^2 * (w(3));

    solut = [a11, a12; ...
             a21, a22] \ [(c^2 / h^2 * sum(I1w * Dxxw * w) - c^2 / h^2 * sum(I1u * Dxxu * u)); ...
                (c^2 / h^2 * sum(I2w * Dxxw * w) - c^2 / h^2 * sum(I2u * Dxxu * u))];
    uNext = uNext + k^2 * (J1u * solut(1) + J2u * solut(2));
    wNext = wNext - k^2 * (J1w * solut(1) + J2w * solut(2));

    scalingU = ones(length(u),1);
    scalingU(end-1:end) = 0.5;
    kinEnergyU(n) = 1/2 * h * sum (scalingU .* (1/k * (u - uPrev)).^2);
    potEnergyU(n) = c^2/(2 * h) * sum((u(2:end-1) - u(1:end-2)) .* (uPrev(2:end-1) - uPrev(1:end-2)));
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0)); % left boundary
    potEnergyU(n) = potEnergyU(n) + c^2/(4 * h) * sum((u(end) - u(end-1)) .* (uPrev(end) - uPrev(end-1))); % right boundary

    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    scalingW = ones(length(w),1);
    scalingW(1:2) = 0.5;

    kinEnergyW(n) = 1/2 * h * sum(scalingW .* (1/k * (w - wPrev)).^2);
    potEnergyW(n) = c^2/(2 * h) * sum((w(3:end) - w(2:end-1)) .* (wPrev(3:end) - wPrev(2:end-1)));
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % left boundary
    potEnergyW(n) = potEnergyW(n) + c^2/(4 * h) * sum((w(2) - w(1)) .* (wPrev(2) - wPrev(1))); % right boundary

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);
        
    totEnergy(n) = totEnergyU(n) + totEnergyW(n);% + connEnergy(n);
    if mod(n, drawSpeed) == 0

        subplot(2,1,1)
        hold off;
        plot(hLocsLeft * N, u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 
%         plot([hLocsLeft, hLocsLeft(end) + h] * N, [0;0;0;0;0;0;-5.55111512312578e-17;-0.146446609406726;-0.353553390593274;-0.500000000000000;-0.500000000000000;-0.353553390593274;-0.146446609406726;0;0;0;0], 'LineWidth' , 2, 'Marker', '.', 'MarkerSize', 20)
        hold on;
        wOffset = 0.00;
        plot(hLocsRight * N, w + wOffset, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
%         plot(z, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 7, 'Color', 'g')
%         plot(hLocsRight * N, [0;5.55111512312578e-1   7;0;0;0.146446609406726;0.353553390593274;0.500000000000000;0.500000000000000;0.353553390593274;0.146446609406726;0;0;0;0;0;0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)
        ylim([-0.6, 0.6])
        grid on;
% %         hold on;
% %         scatter([hLocsLeft(end) + h, hLocsRight(1) - h], [uRightPoint, uLeftPoint]);
% %         scatter((hLocsLeft(end)+h) * N, 0, 'b', 'Marker', 'o', 'Linewidth', 2)
% %         scatter((hLocsRight(1)-h) * N, 0, 'r', 'filled')
% %         legend('$u$', '$w$', '$u_{M+1}$', '$w_{-1}$', 'interpreter', 'latex', 'Fontsize', 16)
%         textOffsetX = 0.1;
%         xlim([12, 19])
%         ylim([-0.1, 0.1] + wOffset * 0.5)
%         textOffsetY = 0.007;
%         text(14 + textOffsetX, -textOffsetY, '$u_{M-2}^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(15 + textOffsetX, -textOffsetY, '$u_{M-1}^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(16 + textOffsetX, -textOffsetY, '$u_M^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(15 + textOffsetX, wOffset + textOffsetY, '$w_0^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(16 + textOffsetX, wOffset + textOffsetY, '$w_1^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(17 + textOffsetX, wOffset + textOffsetY, '$w_2^n$', 'interpreter', 'latex', 'Fontsize', 18)
%         
%         plot([15, 15], [wOffset, 0], '--', 'Color', 'k', 'Linewidth', 2);
%         plot([16, 16], [wOffset, 0], '--', 'Color', 'k', 'Linewidth', 2);
% 
%         text(15 + textOffsetX, wOffset*0.5, '$F_1$', 'interpreter', 'latex', 'Fontsize', 18)
%         text(16 + textOffsetX, wOffset*0.5, '$F_2$', 'interpreter', 'latex', 'Fontsize', 18)
% %         yticklabels([])
%         yticks([])

        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 16
            disp("wait")
        end
        
        subplot(2,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
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
        pause(0.2)
        drawnow;
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
    zPrev = z;
    z = zNext;
    
    outFree(n) = wNext(end - 10);
end
subplot(2,1,1)
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])

c/2