clear all;
close all;

drawSpeed = 100;
fs = 44100;
k = 1 / fs;
lengthSound = fs * 5;
Tinit = 100;
T = Tinit
rho = 7850;
r = 0.0005;
A = pi * r^2;
c = sqrt(Tinit / (rho * A));

h = c*k;
hPrev = h;
N = floor(1/h);
Nprev = N;
h = 1/N;
    
uNext = zeros(N, 1);
u = zeros(N, 1);
u(1:N) = hann(N);
range = 3:N-2;

% u(range) = rand(length(range),1);
uPrev = u;

wait = 1000;
energy = zeros(lengthSound,1);
for n = 1:lengthSound
    if n > wait
        T = Tinit * (1 - (n-wait) / (2 * lengthSound));
    end
    c = sqrt(T/(rho*A));
    
    h = c*k;
%     N = floor(((lengthSound - n)/lengthSound)/h);
    N = floor(1/h);
%     h = 1/N;
    lambdaSq = c^2*k^2/h^2;
    if Nprev ~= N % interpolate
        if N > Nprev
            u = interp1((1:Nprev)*hPrev,u,(1:N)*h,'cubic', 'extrap')';
            uPrev = interp1((1:Nprev)*hPrev,uPrev,(1:N)*h,'cu', 'extrap')';
        else
            u = interp1((1:Nprev)*hPrev,u,(1:N)*h, 'spline', 'extrap')';
            uPrev = interp1((1:Nprev)*hPrev,uPrev,(1:N)*h,'spline', 'extrap')';
        end
        uNext = zeros(N, 1);
    
        
%         uTmp = zeros(N,1);
%         uPrevTmp = zeros(N,1);
%         alphSpace = Nprev / N;
%         % Add the division h1/h2 division here somewhere from stefans book
%         alph = alphSpace;
%         for i = 1:N
%             uTmp(i) = interp (u, i, alph);
%             uPrevTmp(i) = interp (uPrev, i, alph);
%             alph = alph + alphSpace;
%             if alph >= 1
%                 alph = alph - 1;
%             end
%         end
%         uTmp(N) = 0;
%         uPrev(N) = 0;
%         u = uTmp;
%         uPrev = uPrevTmp;
        n
    end
    range = 3:N-2;

    uNext(range) = 2 * u(range) - uPrev(range) + lambdaSq * (u(range-1) - 2 * u(range) + u(range+1));
    if sum(isnan(uNext)) > 0
        disp("wait");
    end
    
    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / 2 * 1/h * sum((u(3:N) - u(2:N-1)) .* (uPrev(3:N) - uPrev(2:N-1)));
    energy(n) = kinEnergy(n) + potEnergy(n);
    
        
    uPrev = u;
    u = uNext;
    out(n) = u(floor(N*2/3));
%     f = test(N);
%     fplot(f)
%     for j = 1:N-1
%         f = u(N-j) + (x - (N-j)*h) * f;
%         fplot(f);
%         drawnow;
%     end
    if mod(n,drawSpeed) == 0
        subplot(2,1,1)
        plot(u);
        ylim([-1, 1])
        subplot(2,1,2)
%         plot(kinEnergy(1:n))
%         hold on;
%         plot(potEnergy(1:n))
%         hold off;
%         plot((energy(8000:n) / energy(8000)) -1)
        plot(energy(1:n))
%         hold off;
    %     fplot(f);
        drawnow;
    end
    Nprev = N;
    hPrev = h;
end

plot(out);

function [ul] = interp (u, n, alph)
    ul = (1-alph) * u(n) + alph * u(n+1);
end