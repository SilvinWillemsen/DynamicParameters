clear all;
close all;
clc;

drawSpeed = 1;
fs = 44100;
k = 1/fs;
lengthSound = fs;

Ninit = 30.0;
N = Ninit;
L = 1;
h = 1/N;

cInit = h/k;
c = cInit;

h = c*k;
N = floor(1/h);
% h = 1/N;
lambdaSq = (c*k/h)^2

overlap = 2;

alf = Ninit - N;

%% initialise states
uNext = zeros(floor(N/2) + ceil((overlap + 1) / 2) - 1, 1);
u = zeros(floor(N/2) + ceil((overlap + 1) / 2) - 1, 1);

if mod(N,2) == 1
    uNext = [uNext; 0];
    u = [u; 0];
end

u(floor(N/5)-4:floor(N/5)+4) = hann (9);
uPrev = u;

wNext = zeros(floor(N/2) + floor((overlap + 1) / 2) - 1, 1);
w = zeros(floor(N/2) + floor((overlap + 1) / 2) - 1, 1);
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

Dxxu(end-1,end-2) = 2 - alf;
Dxxu(end,end-1) = 1 + alf;

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
Dxxw(2,3) = 2 - alf;
Dxxw(1,2) = 1 + alf;

ez = ones(length(z), 1);
Dxxz = spdiags([ez -2*ez ez], -1:1, length(z),length(z));

interpol = "linear";
outFree = zeros(lengthSound, 1);

%% connection  
cub = [alf * (alf - 1) * (alf - 2) / -6, ...
            (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
            alf * (alf + 1) * (alf - 2) / -2, ...
            alf * (alf + 1) * (alf - 1) / 6];
        
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

for n = 1:lengthSound  
    
    hLocsLeft = (1:(length(u))) * h;
    hLocsRight = flip(1 - ((1:(length(w))) * h));
        
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    
    %% orig string
    for i = 1:2
        zNext = 2 * z + lambdaSq * Dxxz * z - zPrev;
        zPrev = z;
        z = zNext;
    end
    a11 = I1u * J1u + I1w * J1w;
    a12 = I1u * J2u + I1w * J2w;
    a21 = I2u * J1u + I2w * J1w;
    a22 = I2u * J2u + I2w * J2w;
    
    solut = [a11, a12; ...
             a21, a22] \ [(c^2 / h^2 * sum(I1w * Dxxw * w) - c^2 / h^2 * sum(I1u * Dxxu * u)); ...
                (c^2 / h^2 * sum(I2w * Dxxw * w) - c^2 / h^2 * sum(I2u * Dxxu * u))];
    uNext = uNext + k^2 * (J1u * solut(1) + J2u * solut(2));
    wNext = wNext - k^2 * (J1w * solut(1) + J2w * solut(2));

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
    if mod(n, drawSpeed) == 0

        subplot(2,1,1)
        hold off;
        plot(hLocsLeft * N, u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 

        hold on;
        wOffset = 0.00;
        plot(hLocsRight * N, w + wOffset, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')

        plot((h/2:h/2:1-h/2) * N, z, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10)

                
        ylim([-0.6, 0.6])
        grid on;

        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 16
            disp("wait")
        end
        
        subplot(2,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
        
        pause(0.2)
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