clear all;
close all;
clc;

drawSpeed = 100000;
drawStart = 1;

fs = 44100;
k = 1/fs;
lengthSound = fs * 2;

Ninit = 35.9;
N = Ninit;  
L = 1;
h = 1/N;

cInit = h/k;
c = cInit;

h = c*k;
N = floor(1/h);
% h = 1/N;
lambdaSq = (c*k/h)^2

alf = Ninit - N;
% alf = 0.5;

%% initialise states
uNext = zeros(floor(N/2), 1);
u = zeros(floor(N/2), 1);

if mod(N,2) == 1
    uNext = [uNext; 0];
    u = [u; 0];
end

u(floor(N/5)-4:floor(N/5)+4) = hann (9);
uPrev = u;

wNext = zeros(floor(N/2), 1);
w = zeros(floor(N/2), 1);
wPrev = w;

zNext = zeros(N*2-1, 1);
z = zeros(N*2-1, 1);
% z(floor(N*2/5)-8:2:floor(N*2/5)+8) = u(floor(N/5)-4:floor(N/5)+4);
% z(floor(N*2/5)-7:2:floor(N*2/5)+7) = 0.5 * (u(floor(N/5)-3:floor(N/5)+4) + u(floor(N/5)-4:floor(N/5)+3));
z(floor(N*2/5)-8:floor(N*2/5)+8) = hann (17);
zPrev = z;

origHLocs = 0:h:1;

eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

% Dxxu(end-1,end-2) = 2 - alf;
Dxxu(end,end-1) = 1;

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
% Dxxw(2,3) = 2 - alf;
Dxxw(1,2) = 1;

ez = ones(length(z), 1);
Dxxz = spdiags([ez -2*ez ez], -1:1, length(z),length(z));

interpol = "cubic";
outFree = zeros(lengthSound, 1);

changeC = false;
for n = 1:lengthSound  
    NPrev = N;
    if changeC
        c = cInit * (1-0.5 * sin(4 * pi * n/lengthSound));
    else
        c = c;
    end
    h = c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambdaSq = c^2 * k^2 / h^2;

    alf = (Ninit - N);
  
    if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
        alf = (1-alf);
        ipRev = [alf * (alf - 1) * (alf - 2) / -6, ...
                    (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                    alf * (alf + 1) * (alf - 2) / -2, ...
                    alf * (alf + 1) * (alf - 1) / 6];

        alf = (1-alf);
    else
        ip = [0, (1-alf), alf, 0];
        alf = (1-alf);
        ipRev = [0, (1-alf), alf, 0];
        alf = (1-alf);
    end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        if mod(N,2) == 1
%             if interpol == "linear"
            uNext = [uNext; (ip(4) * uNext(end-1) + ip(3) * uNext(end) + ip(2) * wNext(1) + ip(1) * wNext(2))];
            u = [u; (ip(4) * u(end-1) + ip(3) * u(end) + ip(2) * w(1) + ip(1) * w(2))];
            uPrev = [uPrev; (ip(4) * uPrev(end-1) + ip(3) * uPrev(end) + ip(2) * wPrev(1) + ip(1) * wPrev(2))];
%             else
%                 uNext = [uNext; ((1-alf) * wNext(1) + alf * wNext(2) - alf * uNext(end)) / (1-alf)];
%                 u = [u; ((1-alf) * w(1) + alf * w(2) - alf * u(end)) / (1-alf)];
%                 uPrev = [uPrev; ((1-alf) * wPrev(1) + alf * wPrev(2) - alf * uPrev(end)) / (1-alf)];
%             
%             end
%             uNext = [uNext; wNext(1)];
%             u = [u; w(1)];
%             uPrev = [uPrev; wPrev(1)];
        else 
            wNext = [(ip(1) * uNext(end-1) + ip(2) * uNext(end) + ip(3) * wNext(1) + ip(4) * wNext(2)); wNext];
            w = [(ip(1) * u(end-1) + ip(2) * u(end) + ip(3) * w(1) + ip(4) * w(2)); w];
            wPrev = [(ip(1) * uPrev(end-1) + ip(2) * uPrev(end) + ip(3) * wPrev(1) + ip(4) * wPrev(2)); wPrev];

%             wNext = [uNext(end); wNext];
%             w = [u(end); w];
%             wPrev = [uPrev(end); wPrev];
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
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
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
    end

    hLocsLeft = (1:(length(u))) * h;
    hLocsRight = flip(1 - ((1:(length(w))) * h));
        
    solut = [1, -ipRev(1); -ip(4), 1] \ [ipRev(2) * w(1) + ipRev(3) * w(2) + ipRev(4) * w(3);
                                        ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
    
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    uNext(end) = uNext(end) + lambdaSq * solut(1);
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + lambdaSq * solut(2);
    
    %% orig string
%     for i = 1:2
%         zNext = 2 * z + lambdaSq * Dxxz * z - zPrev;
%         zPrev = z;
%         z = zNext;
%     end
%     a11 = I1u * J1u + I1w * J1w;
%     a12 = I1u * J2u + I1w * J2w;
%     a21 = I2u * J1u + I2w * J1w;
%     a22 = I2u * J2u + I2w * J2w;
%     
%     solut = [a11, a12; ...
%              a21, a22] \ [(c^2 / h^2 * sum(I1w * Dxxw * w) - c^2 / h^2 * sum(I1u * Dxxu * u)); ...
%                 (c^2 / h^2 * sum(I2w * Dxxw * w) - c^2 / h^2 * sum(I2u * Dxxu * u))];
%     uNext = uNext + k^2 * (J1u * solut(1) + J2u * solut(2));
%     wNext = wNext - k^2 * (J1w * solut(1) + J2w * solut(2));

    %% energies
    
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
    
    outFree(n) = w(end - 5);

    if mod(n, drawSpeed) == 0 && n > drawStart

        subplot(3,1,1)
        hold off;
        plot(hLocsLeft * N, u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 

        hold on;
        wOffset = 0.00;
        plot(hLocsRight * N, w + wOffset, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')

%         plot((h/2:h/2:1-h/2) * N, z, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)

                
        ylim([-0.6, 0.6])
        grid on;

        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 16
            disp("wait")
        end
        
        subplot(3,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
        
        subplot(3,1,3)

        window = 2;
        if n>window
            outfft = fft(outFree(n-window+1 : n));
            data = abs(outfft)  
            semilogy([0:window-1]'*fs/window, data, 'r');
            title("Peak should be at " + num2str(c/2))
            xlim([0 3*cInit/L])
        end
%         pause(0.2)
        drawnow;
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
end

subplot(2,1,1)
hold off
plot(outFree)

subplot(2,1,2)
outfft = fft(outFree);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])

c/2