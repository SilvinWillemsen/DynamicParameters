clear all;
close all;
clc;

drawSpeed = 1;
fs = 44100;
k = 1/fs;

lengthSound = fs;

Ninit = 31.5;
N = Ninit;
if N - floor(N) == 0.5
    virtualFlag = true;
else
    virtualFlag = false;
end

h = 1/N;
c = h/k;

rho = 7850;

r = 0.0005;
A = pi * r^2;
Tinit = c^2 * rho * A;

h = c*k;
N = floor(1/h);
h = 1/N;

lambdaSq = (c*k/h)^2

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
loc = 3/4;
hannRange = floor(N*loc - N/12):floor(N*loc + N/12);
% u(floor(N/2)) = 1;
u(hannRange) = hann(length(hannRange));
% u(floor(N/2)) = 1;
uPrev = u;

origHLocs = 0:h:1-h;

flag = false;
changeT = false;
T = Tinit;

totEnergy1 = [];
kinEnergySave = zeros(N+1, 100);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
potEnergyBound = zeros(lengthSound, 1);

totEnergy = zeros(lengthSound, 1);

interpolationType = "linear";

for n = 1:lengthSound
    NPrev = N;
    if changeT
        T = Tinit * (1-0.5*n/lengthSound);
    else
        T = T;
    end
    
    c = sqrt(T / (rho * A));
    h = c * k;
    
    N = floor(1/h);
    
    hSave(n) = h;
    hLocs = 1:-h:0;
    
    lambdaSq = c^2 * k^2 / h^2;
    range = 2:N;

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        uNext = [0;uNext];
        u = [0;u];
        uPrev = [0;uPrev];
    end
    
    if N < NPrev
        uNext = uNext(2:end);
        u = u(2:end);
        uPrev = uPrev(2:end);
    end
      
    
%     uVirtual = -u(1);
    %% full string
    scaling = ones(N+1, 1);
    scaling(1) = 0.5;
    
    uNext(range) = (2-2*lambdaSq) * u(range) + lambdaSq * (u(range+1) + u(range-1)) - uPrev(range);
    uNext(1) = 2 * u(1) - uPrev(1) + lambdaSq * (u(2) - u(1)) - (1-hLocs(end)/h) * lambdaSq * (u(1) - 0);
    
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2 * h) * sum((u(2:N+1) - u(1:N)) .* (uPrev(2:N+1) - uPrev(1:N)));
    potEnergyBound(n) = T / (2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0));
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + potEnergyBound(n);

    
    
%     totEnergy(n) = kinEnergy(n) + kinEnergyBound(n) + potEnergy(n) + potEnergyBound(n);
    max(kinEnergy(1:n) + potEnergy(1:n) - (kinEnergy(1) + potEnergy(1))) / max(-potEnergyBound(1:n));
    
    if mod(n, drawSpeed) == 0

        subplot(3,1,1)
%         hold off;
        hold off;
        plot(flip(hLocs), u, '-', 'Linewidth', 1, 'Marker', '.', 'MarkerSize', 20, 'Color', [0,0,1])
        hold on; 
        grid on;
        xlim([-2*h, 1])
%         ylim([-1,1])
        set(gca, 'Linewidth', 1)
        
        if n == 13
            disp("wait")
        end
%         scatter([hLocs(end)-hSave(n), hSave(n)-hLocs(end)], [uVirtual, -uVirtual], 400, 'k', 'o')
%         scatter(hLocs(end)-hSave(n), uVirtual, 'k', 'o')
        
        scatter(0, 0, 180, 'k', 'x')
        text(hLocs(end) + h/8, u(1) - 0.05, '$u_0$', 'interpreter', 'latex', 'Fontsize', 26)
%         text(hLocs(end) + h/8 - h, uVirtual-0.05, '$u_{-1}$', 'interpreter', 'latex', 'Fontsize', 26)
        text(hLocs(end) + h + h/8, u(2)-0.05, '$u_1$', 'interpreter', 'latex', 'Fontsize', 26)

%         xlim([-1.5*h, 4*h])
%         plot([hLocs(end) - hSave(n), 0], [uVirtual, 0], '--', 'Linewidth', 2, 'Color', 'b');

        if n == 13
            disp("wait")
        end
        
%         plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual], '--', 'Linewidth', 2, 'Color', 'r');
      
        if n == 13
            disp("wait")
        end
        pause(0.2)
        
        
%         hold on;
%         plot([flip(h1Locs)], u1Plot);
%         scatter(origHLocs, zeros(length(origHLocs),1));
%         scatter(h-hLocs(end), -uVirtual, 18, 'k', 'o')
%         plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual])
% %         xlim([-1.5*h, 4*h])
%         pause (0.25)
        subplot(3,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
%         hold off;
%         plot(totROCEnergy(1:n))
%         hold on;

%         plot(rOCkinEnergy(1:n))
%         plot(rOCpotEnergy(1:n))
%         plot(boundaryEnergy(1:n));
%         subplot(3,1,3)
%         plot(rOCkinEnergy(1:n) - rOCpotEnergy(1:n))
        drawnow;
    end
    
    uPrev = u;
    u = uNext;
    
    out(n) = uNext(floor(N - 10));
end
subplot(2,1,1)
plot(out)
