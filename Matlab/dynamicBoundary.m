clear all;
close all;
clc;

drawSpeed = 100000;
fs = 44100;
k = 1/fs;

lengthSound = fs*2;

% calculate wavespeed from number of points
Ninit = 31.0;
N = Ninit;

h = 1/N;
c = h/k;

% set physical parameters
rho = 7850;
r = 0.0005;
A = pi * r^2;
Tinit = c^2 * rho * A;

% recalculate h
h = c*k;
N = floor(1/h);
h = 1/N;

lambdaSq = (c*k/h)^2

% initialise states (including boundaries)
uNext = zeros(N+1, 1);
u = zeros(N+1, 1);

% excite
loc = 3/4;
hannRange = floor(N*loc - N/12):floor(N*loc + N/12);
u(hannRange) = hann(length(hannRange));
uPrev = u;

origHLocs = 0:h:1-h;

flag = false;
T = Tinit;

totEnergy1 = [];
kinEnergySave = zeros(N+1, 100);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
potEnergyBound = zeros(lengthSound, 1);
potEnergyBound2 = zeros(lengthSound, 1);

totEnergy = zeros(lengthSound, 1);

interpolationType = "linear";
uVirtualPrev = 0;
potEnergyRange = 1:N;
changeT = false;
L = 1;
for n = 1:lengthSound
    NPrev = N;
    if changeT
        T = Tinit * (1-0.5*n/lengthSound);
%         T = Tinit * ((1 - (n-10)/10) + 0.5^2 * (n-10)/10);
%         L = 1 + (n-10)/10;
    else
        T = T;
    end
    
    c = sqrt(T / (rho * A));
    h = c * k;
    
    N = round(10000 * L/h) / 10000;
    
    hSave(n) = h;
    hLocs = L:-h:0;
    
    lambdaSq =  c^2 * k^2 / h^2;
    range = 2:N;
  
    if abs(N - NPrev) > 1
        disp('too fast')
    end
    while N > length(uNext) - 1
        uNext = [0;uNext];
        u = [0;u];
        uPrev = [0;uPrev];
    end
    
    if N < NPrev
        uNext = uNext(2:end);
        u = u(2:end);
        uPrev = uPrev(2:end);
    end
      
    if hLocs(end) >= hSave(n)/2
        alpha = (hSave(n) - hLocs(end)) / hLocs(end);
        uVirtual = -alpha * u(1);
    else
        alpha = (2*hLocs(end)) / hSave(n);
        uVirtual = -(alpha * u(1) + (1-alpha) * u(2));
    end
    
    %% full string
    uNext(range) = (2-2*lambdaSq) * u(range) + lambdaSq * (u(range+1) + u(range-1)) - uPrev(range);
    uNext(1) = 2 * u(1) - uPrev(1) + lambdaSq * (u(2) - 2 * u(1) + uVirtual);

%     kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u(range) - uPrev(range))).^2);
%     kinEnergyBound(n) = rho * A / 2 * h * (1/k * (u(1) - uPrev(1)))^2;

    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2 * h) * sum((u(potEnergyRange+1) - u(potEnergyRange)) .* (uPrev(potEnergyRange+1) - uPrev(potEnergyRange)));
%     potEnergyBound(n) = (1+alpha) * T / (2 * h) * (u(1) - 0) .* (uPrev(1) - 0);
    potEnergyBound(n) = T / (2 * h * (Ninit - floor(Ninit))) * (u(1) - 0) * (uPrev(1) - 0);
%     potEnergyBound(n) = T / (2 * h * (1+alpha)) * (u(1) - uVirtual) * (uPrev(1) - uVirtualPrev);
    potEnergyBound2(n) = T / (2 * h) * (0 - uVirtual) * (0 - uVirtualPrev);
%     potEnergyBound(n) = (1 + alpha) * T / (2 * h) * sum((-uVirtual/alpha - 0) .* (-uVirtualPrev/alpha - 0));
    
    rOCkinEnergy(n) = rho * A * h / (2*k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = -T / (2 * k * h) * sum((u(range+1) - 2 * u(range) + u(range-1)) .* (uNext(range) - uPrev(range)));
    rOCboundEnergy(n) = -T / (2 * k * h) * (u(2) - 2 * u(1) + uVirtual) * (uNext(1) - uPrev(1));
    rOCTotEnergy(n) = rOCkinEnergy(n) + rOCpotEnergy(n) + rOCboundEnergy(n);
    
    rOCpotTest(n) = T / (2 * k * h) * sum((uNext(potEnergyRange+1) - uNext(potEnergyRange) - uPrev(potEnergyRange+1) + uPrev(potEnergyRange)) .* (u(potEnergyRange+1) - u(potEnergyRange)));
    rOCpotTestBound(n) =  T / (2 * k * h) * (uNext(1) - uPrev(1)) * (u(1) - uVirtual);

    
    rOCTotEnergyTest(n) = rOCkinEnergy(n) + rOCpotTest(n) + rOCpotTestBound(n);
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + potEnergyBound(n) + potEnergyBound2(n);

%     kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2 * h) * sum((u(2:N+1) - u(1:N)) .* (uPrev(2:N+1) - uPrev(1:N)));
    potEnergyBound(n) = 2 * T / (2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0));
    scaling = ones(N+1, 1);
    scaling(1) = 0.5;
    
%     rOCpotRange = 2:N;
%     rOCkinEnergy(n) = h * rho * A / (2 * k^3) * sum(scaling .* (uNext - 2 * u + uPrev) .* (uNext - uPrev));
%     rOCpotEnergy(n) = h * T / (2*k*h^2) * sum((u(rOCpotRange+1) - 2 * u(rOCpotRange) + u(rOCpotRange-1))...
%         .* (uNext(rOCpotRange) - uPrev(rOCpotRange)));
%     boundaryEnergy(n) = -T / (4 * k * h)  * (uNext(1) - uPrev(1)) * (u(1) + u(2));
%     
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + potEnergyBound(n);
        
    if mod(n, drawSpeed) == 0

%         subplot(3,1,1)
%         hold off;
        hold off;
        plot(flip(hLocs), u, '-', 'Linewidth', 1, 'Marker', '.', 'MarkerSize', 20, 'Color', [0,0,1])
        hold on; 
        grid on;
        xlim([0, L])
        ylim([-.6,.6])
%         title("Time step $n = $" + num2str(n) + "  Wave speed $c = $" + num2str(c) + "  No. points $N = $" + num2str(N), 'interpreter', 'latex');
     
        set(gca, 'Fontsize', 16);
        set(gca, 'Linewidth', 1)
        
        if n == 13
            disp("wait")
        end
%         scatter([hLocs(end)-hSave(n), hSave(n)-hLocs(end)], [uVirtual, -uVirtual], 400, 'k', 'o')
%         scatter(hLocs(end)-hSave(n), uVirtual, 'k', 'o')
%         
%         scatter(0, 0, 180, 'k', 'x')
%         text(hLocs(end) + h/8, u(1) - 0.05, '$u_0$', 'interpreter', 'latex', 'Fontsize', 26)
%         text(hLocs(end) + h/8 - h, uVirtual-0.05, '$u_{-1}$', 'interpreter', 'latex', 'Fontsize', 26)
%         text(hLocs(end) + h + h/8, u(2)-0.05, '$u_1$', 'interpreter', 'latex', 'Fontsize', 26)

%         xlim([-1.5*h, 4*h])
%         plot([hLocs(end) - hSave(n), 0], [uVirtual, 0], '--', 'Linewidth', 2, 'Color', 'b');

        if n == 13
            disp("wait")
        end
        
%         plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual], '--', 'Linewidth', 2, 'Color', 'r');
      
        if n == 13
            disp("wait")
        end
        pause(0.5)
        
        
%         hold on;
%         plot([flip(h1Locs)], u1Plot);
%         scatter(origHLocs, zeros(length(origHLocs),1));
%         scatter(h-hLocs(end), -uVirtual, 18, 'k', 'o')
%         plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual])
% %         xlim([-1.5*h, 4*h])
%         pause (0.25)
%         subplot(3,1,2)
%         plot(totEnergy(1:n) / totEnergy(1) - 1)
%         hold off;
%         plot(totROCEnergy(1:n))
%         hold on;
%         hold off;
%         plot(rOCkinEnergy(1:n))
%         hold on;
%        
%         plot(rOCpotTest(1:n))
%        plot(rOCTotEnergyTest)
      
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
% plot(out)

outfft = fft(out);
plot([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c/L])
c/2