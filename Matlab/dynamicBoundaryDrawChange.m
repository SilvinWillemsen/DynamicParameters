clear all;
close all;
clc;

drawSpeed = 1;
fs = 44100;
k = 1/fs;

lengthSound = fs*2;

Ninit = 30.1;
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
changeT = true;
T = Tinit;

totEnergy1 = [];
kinEnergySave = zeros(N+1, 100);

kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
potEnergyBound = zeros(lengthSound, 1);

totEnergy = zeros(lengthSound, 1);

interpolationType = "linear";
uVirtualPrev = 0;
for n = 1:lengthSound
    NPrev = N;
    if changeT
        T = Tinit;% * (1-0.5*n/lengthSound);
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
      
    if interpolationType == "linear"
        if hLocs(end) >= hSave(n)/2
            alpha = (hSave(n) - hLocs(end)) / hLocs(end);
            uVirtual = -alpha * u(1);
            flag = true;
        else
            alpha = (2*hLocs(end)) / hSave(n);
            uVirtual = -(alpha * u(1) + (1-alpha) * u(2));
        end
    elseif interpolationType == "cubic"
        res = interp1([0, hLocs(end:-1:end-2)],[0, u(1:3)'], -(hLocs(end)-h), 'spline');
        uVirtual = -res;
    end
    
    %% full string
    uNext(range) = (2-2*lambdaSq) * u(range) + lambdaSq * (u(range+1) + u(range-1)) - uPrev(range);
    uNext(1) = 2 * u(1) - uPrev(1) + lambdaSq * (u(2) - 2 * u(1) + uVirtual);

%     kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u(range) - uPrev(range))).^2);
%     kinEnergyBound(n) = rho * A / 2 * h * (1/k * (u(1) - uPrev(1)))^2;

    kinEnergy(n) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = T / (2 * h) * sum((u(2:N+1) - u(1:N)) .* (uPrev(2:N+1) - uPrev(1:N)));
    potEnergyBound(n) = (1 + alpha) * T / (2 * h) * (u(1) - 0) .* (uPrev(1) - 0);
    potEnergyBound(n) = T / (2 * h) * (u(1) - 0) * (uPrev(1) - 0);
    potEnergyBound2(n) = T / (2 * h * alpha) * (0 - uVirtual) * (0 - uVirtualPrev);
%     potEnergyBound(n) = (1 + alpha) * T / (2 * h) * sum((-uVirtual/alpha - 0) .* (-uVirtualPrev/alpha - 0));
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + potEnergyBound(n) + potEnergyBound2(n);

    scaling = ones(N+1, 1);
    scaling(1) = 0.5;
    
%     totEnergy(n) = kinEnergy(n) + kinEnergyBound(n) + potEnergy(n) + potEnergyBound(n);
%     max(kinEnergy(1:n) + potEnergy(1:n) - (kinEnergy(1) + potEnergy(1))) / max(-potEnergyBound(1:n));
    
    if mod(n, drawSpeed) == 0
        stopIdx = 1;
%         subplot(3,1,1)
%         hold off;
    for ii = 1:5
        Ntest = 30.0 + 0.2499999 * (ii-1)

        hSave(n) = 1/Ntest;
        h = hSave(n);
        hLocs = 1:-h:0;
        hold on
%         hold off;
        plot(flip(hLocs), u + ii, '-', 'Linewidth', 1, 'Marker', '.', 'MarkerSize', 20, 'Color', [0,0,1])
        hold on; 
        grid on;
        xlim([0, 1])
%         ylim([-.6,.6])
        set(gca, 'Linewidth', 1)
        
        if n == stopIdx
            disp("wait")
        end
%         scatter([hLocs(end)-hSave(n), hSave(n)-hLocs(end)], [uVirtual, -uVirtual], 400, 'k', 'o')
        
        
        scatter(0, 0 + ii, 180, 'k', 'x')
        xoffset = h/15;
        yoffset = 0.15;
        if ii == 2
            test = 1;
            plot([0, h - hLocs(end)], [ii+0.2, ii+0.2], 'r', 'Linewidth', test)
            plot([0, 0], [ii+0.1, ii+0.3], 'r', 'Linewidth', test)
            plot([h-hLocs(end), h-hLocs(end)], [ii+0.1, ii+0.3], 'r', 'Linewidth', test)
            text(0.5 * (h - hLocs(end)) - xoffset, ii + 0.4, '$h_{\textrm{\fontsize{7}{0}\selectfont I}}$', 'interpreter', 'latex', 'Fontsize', 20, 'Color', 'r')

        end
        if ii ~= 5
            text(hLocs(end) + xoffset, u(1) - yoffset + ii, '$u_0$', 'interpreter', 'latex', 'Fontsize', 26)
            text(hLocs(end) + xoffset - h, uVirtual-yoffset + ii, '$u_{-1}$', 'interpreter', 'latex', 'Fontsize', 26)
            text(hLocs(end) + h + xoffset, u(2)-yoffset + ii, '$u_1$', 'interpreter', 'latex', 'Fontsize', 26)
            text((h-hLocs(end)) + xoffset, uVirtual + yoffset*1.3 + ii, '$u_{\textrm{\fontsize{7}{0}\selectfont I}}$', 'interpreter', 'latex', 'Fontsize', 26, 'Color', 'r')

            plot([hLocs(end) - hSave(n), hLocs(end)], [uVirtual, 0] + ii, '--', 'Linewidth', 2, 'Color', 'b');
            scatter(hLocs(end)-hSave(n), uVirtual + ii, 'k', 'o')
        else
            scatter(hLocs(end)-hSave(n), uVirtual + ii, 400, 'b', '.');
            text(hLocs(end) + xoffset, u(1) - yoffset + ii, '$u_1$', 'interpreter', 'latex', 'Fontsize', 26)
            text(hLocs(end) + xoffset - h, uVirtual-yoffset + ii, '$u_0$', 'interpreter', 'latex', 'Fontsize', 26)
            text(hLocs(end) + h + xoffset, u(2)-yoffset + ii, '$u_2$', 'interpreter', 'latex', 'Fontsize', 26)
            plot([hLocs(end) - hSave(n), hLocs(end)], [uVirtual, 0] + ii, 'Linewidth', 1, 'Color', 'b');

        end   
        xlim([-1.5*h, 2.5*h])
        ylim([0.5, 5.5])
        if n == stopIdx
            disp("wait")
        end
        
        plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual] + ii, '--', 'Linewidth', 2, 'Color', 'r');
      
        if n == stopIdx
            disp("wait")
        end
%         pause(0.2)
        
        
%         hold on;
%         plot([flip(h1Locs)], u1Plot);
%         scatter(origHLocs, zeros(length(origHLocs),1));
        if ii ~= 5
            scatter(h-hLocs(end), -uVirtual + ii, 400, 'r', '.')
        end
%         plot([hLocs(end) - hSave(n), h-hLocs(end)], [uVirtual, -uVirtual])
% % %         xlim([-1.5*h, 4*h])
% %         pause (0.25)
%         subplot(3,1,2)
%         plot(totEnergy(1:n) / totEnergy(1) - 1)
%         hold off;
%         plot(totROCEnergy(1:n))
%         hold on;

%         plot(rOCkinEnergy(1:n))
%         plot(rOCpotEnergy(1:n))
%         plot(boundaryEnergy(1:n));
%         subplot(3,1,3)
%         plot(rOCkinEnergy(1:n) - rOCpotEnergy(1:n))
%         drawnow;
    end
    xticks([0])
    xticklabels({'$u_{\textrm{\fontsize{7}{0}\selectfont B}}$'})
%     gca
    set(gca, 'TickLabelInterpreter', 'latex')
    yticks([1:5])
    yticklabels({'','','','',''})

    ylabel('Time')
    set(gca, 'Fontsize', 16)
    
    end
    uVirtualPrev = uVirtual;
    uPrev = u;
    u = uNext;
    
    out(n) = uNext(floor(N - 10));
end
subplot(2,1,1)
plot(out)
