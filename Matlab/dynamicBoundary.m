clear all;
close all;
clc;

upsample = 4;

drawSpeed = 1;
fs = 44100;
fs1 = fs * upsample;
k = 1/fs;
k1 = 1/fs1;

lengthSound = fs*2;

Ninit = 30.5;
N = Ninit;
if N - floor(N) == 0.5
    virtualFlag = true;
else
    virtualFlag = false;
end

h = 1/N;
h1 = 1/(N*upsample);

c = h/k;

rho = 7850;

r = 0.0005;
A = pi * r^2;
Tinit = c^2 * rho * A;

h = c*k;
h1 = c*k1;

N = floor(1/h);
N1 = floor(1/h1);

h = 1/N;
h1 = 1/N1;

lambdaSq = (c*k/h)^2

uNext = zeros(N+1, 1);
u = zeros(N+1, 1);
loc = 3/4
hannRange = floor(N*loc - N/12):floor(N*loc + N/12);
% u(floor(N/2)) = 1;
u(hannRange) = hann(length(hannRange));
% u(floor(N/2)) = 1;
uPrev = u;

u1Next = zeros(N1+1, 1);
u1 = zeros(N1+1, 1);
hannRange1 = floor(N1*loc - N1/12):floor(N1*loc + N1/12);
u1(hannRange1) = hann(length(hannRange1));
% u1(floor(N1/2)) = 1;

u1Prev = u1;

origHLocs = 0:h:1-h;
origH1Locs = 0:h1:1-h1;

flag = false;
changeT = true;
T = Tinit;

totEnergy1 = [];
kinEnergySave = zeros(N+1, 100);
for n = 1:lengthSound
    NPrev = N;
    if changeT
        T = Tinit;% * (1-0.5*n/lengthSound);
    else
        T = T;
    end
    
    c = sqrt(T / (rho * A));
    h = c * k;
%     h1 = c * k1;
    
    N = floor(1/h);
%     N1 = floor(1/h1);
    
    hSave(n) = h;
    hLocs = 1:-h:0;
    
%     h1Save(n) = h1;
%     h1Locs = 1:-h1:0;

    lambdaSq =  c^2 * k^2 / h^2;
    range = 2:N;
%     range1 = 2:N1;

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        uNext = [0;uNext];
        u = [0;u];
        uPrev = [0;uPrev];
%         u1Next = [0;u1Next];
%         u1 = [0;u1];
%         u1Prev = [0;u1Prev];
    end
    
    if N < NPrev
        uNext = uNext(2:end);
        u = u(2:end);
        uPrev = uPrev(2:end);
%         u1Next = u1Next(2:end);
%         u1 = u1(2:end);
%         u1Prev = u1Prev(2:end);
    end
      
    if hLocs(end) >= hSave(n)/2
        alpha = (hSave(n) - hLocs(end)) / hLocs(end);
        uVirtual = -alpha * u(1);
        flag = true;
    else
        alpha = (2*hLocs(end)) / hSave(n);
        uVirtual = -(alpha * u(1) + (1-alpha) * u(2));
    end
    
    if virtualFlag
        uVirtual = -u(1);
    end
    
    %% full string
    uNext(range) = (2-2*lambdaSq) * u(range) + lambdaSq * (u(range+1) + u(range-1)) - uPrev(range);
    uNext(1) = 2 * u(1) - uPrev(1) + lambdaSq * (u(2) - 2 * u(1) + uVirtual);
%     for m = 1:upsample
%         u1Next(range1) = (2-2*lambdaSq) * u1(range1) + lambdaSq * (u1(range1+1) + u1(range1-1)) - u1Prev(range1);
%  
%         kinEnergy1 = 1 / 2 * h1 * sum((1/k1 * (u1 - u1Prev)).^2);
%         potEnergy1 = c^2 / (2 * h1) * sum((u1(range1) - u1(range1-1)) .* (u1Prev(range1) - u1Prev(range1-1)));
%         totEnergy1 = [totEnergy1; kinEnergy1+potEnergy1];
%         if m == 1
%             u1Plot = u1Next;
%         end
%         u1Prev = u1;
%         u1 = u1Next;
%     end    
%     
    rOCpotRange = 2:N;
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u(range) - uPrev(range))).^2);
    kinEnergyBound(n) = rho * A / 2 * h * (1/k * (u(1) - uPrev(1)))^2;

%     kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2 * h) * sum((u(2:N+1) - u(1:N)) .* (uPrev(2:N+1) - uPrev(1:N)));
    scaling = ones(N+1, 1);
    scaling(1) = 0.5;
    rOCkinEnergy(n) = h * rho * A / (2 * k^3) * sum(scaling .* (uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = h * T / (2*k*h^2) * sum((u(rOCpotRange+1) - 2 * u(rOCpotRange) + u(rOCpotRange-1))...
        .* (uNext(rOCpotRange) - uPrev(rOCpotRange)));
    boundaryEnergy(n) = -T / (4 * k * h)  * (uNext(1) - uPrev(1)) * (u(1) + u(2));
    
    totEnergy(n) = kinEnergy(n) + kinEnergyBound(n) + potEnergy(n);
    
    totROCEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - boundaryEnergy(n);
    
    if mod(n, drawSpeed) == 0

        subplot(3,1,1)
%         hold off;
        hold off;
        plot(flip(hLocs), u, '-', 'Linewidth', 1, 'Marker', '.', 'MarkerSize', 20, 'Color', [0,0,1])
        hold on; 
        grid on;
        xlim([0, 1])
        ylim([-.6,.6])
        set(gca, 'Linewidth', 1)
        
        if n == 13
            disp("wait")
        end
%         scatter([hLocs(end)-hSave(n), hSave(n)-hLocs(end)], [uVirtual, -uVirtual], 400, 'k', 'o')
        scatter(hLocs(end)-hSave(n), uVirtual, 'k', 'o')
        
        scatter(0, 0, 180, 'k', 'x')
        text(hLocs(end) + h/8, u(1) - 0.05, '$u_0$', 'interpreter', 'latex', 'Fontsize', 26)
        text(hLocs(end) + h/8 - h, uVirtual-0.05, '$u_{-1}$', 'interpreter', 'latex', 'Fontsize', 26)
        text(hLocs(end) + h + h/8, u(2)-0.05, '$u_1$', 'interpreter', 'latex', 'Fontsize', 26)

%         xlim([-1.5*h, 4*h])
        plot([hLocs(end) - hSave(n), 0], [uVirtual, 0], '--', 'Linewidth', 2, 'Color', 'b');

        if n == 13
            disp("wait")
        end
        
        plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual], '--', 'Linewidth', 2, 'Color', 'r');
      
        if n == 13
            disp("wait")
        end
        pause(0.1)
        
        
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
