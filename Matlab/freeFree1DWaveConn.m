clear all;
close all;
clc;

drawSpeed = 10;
fs = 44100;
k = 1/fs;
lengthSound = fs*2;

N = 50.25;

h = 1/N;

c = h/k;
rho = 7850;
Einit = c^2 * rho;

h = c*k;
N = floor(1/h);
h = 1/N;
lambdaSq = (c*k/h)^2

u1Next = zeros(N+1, 1);
u1 = zeros(N+1, 1);
u1(floor(N/4)-2:floor(N/4)+2) = hann(5);
u1Prev = u1;

u2Next = zeros(floor(N/2) + 1, 1);
u2 = zeros(floor(N/2) + 1, 1);
u2(floor(N/4)-2:floor(N/4)+2) = hann(5);
u2Prev = u2;

u3Next = zeros(floor(N/2) + 1, 1);
u3 = zeros(floor(N/2) + 1, 1);
% u3(floor(N/4)-5:floor(N/4)+5) = hann(11);
u3Prev = u3;

uNext = zeros(N, 1);
u = zeros(N, 1);
u(floor(N/4)-2:floor(N/4)+2) = hann(5);
uPrev = u;

halfRange = 2:length(u2)-1;

origHLocs = 0:h:1;

flag = false;
changeE = true;
E = Einit;
for n = 1:lengthSound
    NPrev = N;
    if changeE
        E = Einit * (1-0.5*n/lengthSound);
    else
        E = E;
    end
    K = E / h;
    M = rho * h;
    h = k * sqrt(E / rho);
    N = floor(1/h);
    hSave(n) = h;
    hLocs = 1:-h:0;

    lambdaSq = E / rho * k^2 / h^2;
    range = 2:N;

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    if N > NPrev
        uNext = [0;uNext];
        u = [0;u];
        uPrev = [0;uPrev];
        u1Next = [0;u1Next];
        u1 = [0;u1];
        u1Prev = [0;u1Prev];
    end
    
    if N < NPrev
        uNext = uNext(2:end);
        u = u(2:end);
        uPrev = uPrev(2:end);
        u1Next = u1Next(2:end);
        u1 = u1(2:end);
        u1Prev = u1Prev(2:end);
    end
%     %% mass-spring string
%     for l = 2:N-1
%         uNext(l) = 2 * u(l) - uPrev(l) + k^2 * K/M * (u(l+1) - u(l)) - k^2 * K/M * (u(l) - u(l-1));
%     end
%         
    if hLocs(end) >= hSave(n)/2
        alpha = (hSave(n) - hLocs(end)) / hLocs(end);
        u1Virtual = -alpha * u1(1); %0 = boundary;
        flag = true;
    else
        alpha = (2*hLocs(end)) / hSave(n);
        u1Virtual = -(alpha * u1(1) + (1-alpha) * u1(2));
    end
%         if flag && floor(alpha * 1000) == 0
%             changeE = false;
%         end

%     u1Virtual = -tand(atand(u1(1)/(alpha*h)))*(1-alpha*h);
%     u1Virtual = u1(2);
%     if u1Virtual ~= 0
%         disp("wait")
%     end
    %% full string
    u1Next(range) = (2-2*lambdaSq) * u1(range) + lambdaSq * (u1(range+1) + u1(range-1)) - u1Prev(range);
%     u1Next(1) = (2 * u1(1) - u1Prev(1) + 2 * lambdaSq * (u1(2) - u1(1)) + alpha * omega0^2 * k^2 / (2 * h) * uPrev(1) - (1-alpha) * lambdaSq * u1(2)) / (1 + alpha * omega0^2 * k^2 / (2 * h));
    u1Next(1) = 2 * u1(1) - u1Prev(1) + lambdaSq * (u1(2) - 2 * u1(1) + u1Virtual);
    
    kinEnergy(n) = 1 / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy(n) = c^2 / (2 * h) * sum((u1(2:N+1) - u1(1:N)) .* (u1Prev(2:N+1) - u1Prev(1:N)));
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n);
%     %% left half string
%     u2Next(halfRange) = (2-2*lambdaSq) * u2(halfRange) + lambdaSq * (u2(halfRange+1) + u2(halfRange-1)) - u2Prev(halfRange);
%     u2Next(end) = (2-2*lambdaSq) * u2(end) + lambdaSq * (2 * u2(end-1)) - u2Prev(end);
%     
%     %% right half string
%     u3Next(halfRange) = (2-2*lambdaSq) * u3(halfRange) + lambdaSq * (u3(halfRange+1) + u3(halfRange-1)) - u3Prev(halfRange);
%     u3Next(1) = (2-2*lambdaSq) * u3(1) + lambdaSq * (2 * u3(2)) - u3Prev(1);
%         
%     %% connection
%     F = h * (u3Next(1) - u2Next(end)) / 2;
%     
%     u2Next(end) = u2Next(end) + 1/h * F;
%     u3Next(1) = u3Next(1) - 1/h * F;
    
    if mod(n, drawSpeed) == 0
%         hold off;
%         scatter(flip(hLocs), zeros(length(hLocs),1));
%         hold on;
%         scatter(origHLocs, zeros(length(origHLocs),1));
        subplot(2,1,1)
%         plot(uNext, 'o-');
%         subplot(2,1,2)
        hold off;
        plot([hLocs(end)-hSave(n), flip(hLocs)], [u1Virtual; u1], 'o-')
        ylim([-1,1])
        hold on;
        scatter(origHLocs, zeros(length(origHLocs),1));
        scatter(0, 0, 180, 'k', 'x')
        scatter(h-hLocs(end), -u1Virtual, 18, 'k', 'o')
        plot([hLocs(end) - hSave(n), h-hLocs(end)], [u1Virtual, -u1Virtual])
%         subplot(3,1,3)
%         plot(flip(hLocs), [u2Next;u3Next], 'o-')
    xlim([-1.5*h, 4*h])
    pause (0.25)
        subplot(2,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
        drawnow;
    end
    uPrev = u;
    u = uNext;
    
    u1Prev = u1;
    u1 = u1Next;
    
    u2Prev = u2;
    u2 = u2Next;
    
    u3Prev = u3;
    u3 = u3Next;
    out1(n) = u1Next(floor(N - 10));
    out(n) = uNext(floor(N - 10));
%     sum(u1-[u2;u3(2:end)])
end
subplot(2,1,1)
plot(out)
subplot(2,1,2)
plot(out1)