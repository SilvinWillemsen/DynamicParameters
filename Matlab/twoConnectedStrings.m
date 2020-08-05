clear all;
close all;

fs = 44100;
lengthSound = fs;
k = 1/fs;

N = 30;
h = 1/N;
c = h/k;

uNext = zeros(N+1, 1);
u = uNext;
u(N/4-4:N/4+4) = hann(9)
uPrev = u;

wNext = zeros(N+1, 1);
w = wNext;
% w(N/4-4:N/4+4) = hann(9)
wPrev = w;
lambdaSq = (c * k / h)^2;

range = 2:N;
energyRange = 1:N;

alpha = 0.01;

Iu = zeros(1, N + 1);
Iu(floor(N/2)) = (1-alpha);
Iu(floor(N/2) + 1) = alpha;

Iw = zeros(1, N + 1);
Iw(floor(N/2)) = (1-alpha);
Iw(floor(N/2) + 1) = alpha;

omega0 = 10000;
for n = 1:lengthSound
    
        %calculate schemes
        uNext(range) = 2 * u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2 * u(range) + u(range-1));
        wNext(range) = 2 * w(range) - wPrev(range) + lambdaSq * (w(range+1) - 2 * w(range) + w(range-1));
        
        % calculate forces
        eta = Iu * u - Iw * w;
        etaPrev = Iu * uPrev - Iw * wPrev;
        F = (Iw * wNext - Iu * uNext - etaPrev) / (2 * k^2 / h + 2 / omega0^2);
        
        % add forces
        uNext = uNext + k^2/h * Iu' * F;
        wNext = wNext - k^2/h * Iu' * F;

        
        % energies
        Hu(n) = 1/2 * h * sum((1/k * (u - uPrev)).^2) ... % kinEnergy
            + c^2 / (2 * h) * sum((u(energyRange + 1) - u(energyRange)) .* (uPrev(energyRange + 1) - uPrev(energyRange))); % potEnergy
        Hw(n) = 1/2 * h * sum((1/k * (w - wPrev)).^2) ... % kinEnergy
            + c^2 / (2 * h) * sum((w(energyRange + 1) - w(energyRange)) .* (wPrev(energyRange + 1) - wPrev(energyRange))); % potEnergy

        conn(n) = omega0^2 / 2 * 1/2 * (eta^2 + etaPrev^2);
        
        totEnergy(n) = Hu(n) + Hw(n);% + conn(n);
        
        subplot(2,1,1)
        hold off
        plot(u, '-o');
        hold on;
        plot(w, '-o');
        
        subplot(2,1,2)
        hold off;
        plot(-(totEnergy(1:n) - totEnergy(1)));
        hold on;
        plot(conn(1:n))
        pause(0.2)
        drawnow;
        
        uPrev = u;
        u = uNext;
        
        wPrev = w;
        w = wNext;
end