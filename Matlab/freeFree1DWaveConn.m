clear all;
close all;
clc;

drawSpeed = 1;
drawStart = 0;
drawThings = true;

fs = 44100;             % Sample rate
k = 1/fs;               % Time step
lengthSound = fs * 1;       % Length of the simulation

startSample = 0;

Ninit = 30 * fs / 44100;           % edit how many points you want
h = 1/Ninit;
Nend = 31 * fs / 44100;
cEnd = 1/(Nend*k);
cInit = h/k;            % calculate wave speed
c = cInit;

h = c*k;                % calculate h from wavespeed
N = floor(1/h);         % calculate points from h

lambdaSq = (c*k/h)^2    % should always be 1 as h is not recalculated

alf = Ninit - N;        % fractional remainder for the grid point
% alf = 0;
%% initialise states
uNext = zeros(ceil(N/2), 1);
u = zeros(ceil(N/2), 1);

excitationWidth = 0.2;
excitationLoc = 0.2;
loc = excitationLoc * N;
width = max (4.0, excitationWidth * N);
raisedCosStart = floor (loc - width * 0.5);
raisedCosEnd = floor (loc + width * 0.5);

u(raisedCosStart : raisedCosEnd) = 0.5 * (1 - cos (2.0 * pi * ((raisedCosStart:raisedCosEnd)-raisedCosStart) / width));
% u(floor(N/5)-4:floor(N/5)+4) = hann(9); % use hann window for excitation
% u = rand(length(u), 1);
uPrev = u;

wNext = zeros(floor(N/2), 1);
w = zeros(floor(N/2), 1);
% w(floor(4*/5)-4:floor(4*N/5)+4) = hann(9); % use hann window for excitation

% w = rand(length(w), 1);
% w(1) = u(end);
wPrev = w;

% initialise laplacians
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

interpol = "cubic";
outFree = zeros(lengthSound, 1);

changeC = true; % set to true for dynamic changes in wavespeed

%% recording
M(500) = struct('cdata',[],'colormap',[]);
frame = 1;
filmFlag = true;
interpolatedPoints = [0; 0];

%% plotting

cVec = linspace(cInit, cEnd, lengthSound);
figure('Position', [100, 100, 500, 210])

uSave = [];
wSave = [];
MuSave = [];
MwSave = [];

NVec = linspace(Ninit, Nend, lengthSound);

for n = 1:lengthSound  
 
%     % change wave speed
    if changeC
        c = cVec(n);
    elseif changeN
        Ninit = NVec(n);
        h = 1/Ninit;
        c = h/k;
    else
        c = c;
    end
%     cSave(n) = c;
    
    % save previous state for comparison later
    NPrev = N;

  
    cSave(n) = c;
    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambdaSq = c^2 * k^2 / h^2;

    alf = (Ninit - N);
  alfSave(n) = alf;
    % calculate interpolator
    if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
    else
        ip = [0, (1-alf), alf, 0];
    end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    
    % add point if N^n > N^{n-1}
    if N > NPrev
        if mod(N,2) == 1
            uNext = [uNext; (ip(4) * uNext(end-1) + ip(3) * uNext(end) + ip(2) * wNext(1) + ip(1) * wNext(2))];
            u = [u; (ip(4) * u(end-1) + ip(3) * u(end) + ip(2) * w(1) + ip(1) * w(2))];
            uPrev = [uPrev; (ip(4) * uPrev(end-1) + ip(3) * uPrev(end) + ip(2) * wPrev(1) + ip(1) * wPrev(2))];
        else 
            wNext = [(ip(1) * uNext(end-1) + ip(2) * uNext(end) + ip(3) * wNext(1) + ip(4) * wNext(2)); wNext];
            w = [(ip(1) * u(end-1) + ip(2) * u(end) + ip(3) * w(1) + ip(4) * w(2)); w];
            wPrev = [(ip(1) * uPrev(end-1) + ip(2) * uPrev(end) + ip(3) * wPrev(1) + ip(4) * wPrev(2)); wPrev];
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));
    end
    
    % remove point if N^n < N^{n-1}
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
     
    % calculate interpolated points
    interpolatedPointsPrev = interpolatedPoints;
    interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                        ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
    
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    uNext(end) = uNext(end) + lambdaSq * interpolatedPoints(1);
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + lambdaSq * interpolatedPoints(2);

    if n > startSample && n < startSample + 100
        uSave = [uSave; u];
        wSave = [wSave; w];
        MuSave = [MuSave; length(u)];
        MwSave = [MwSave; length(w)];
    end
    %% energies
    
    scalingU = ones(length(u),1);
    scalingU(end) = 0.5 * (1 + alf);
%     kinEnergyU(n) = 1/2 * h * sum (scalingU .* (1/k * (u - uPrev)).^2);
    kinEnergyU(n) = 1/2 * h * sum ((1/k * (u(1:end-1) - uPrev(1:end-1))).^2);
    connKinEnergyU(n) = 1/2 * h * scalingU(end) * (1/k * (u(end) - uPrev(end))).^2;
    
    potEnergyU(n) = c^2/(2 * h) * sum((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)));
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0)); % left boundary
%     potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((w(2) - u(end)) .* (wPrev(2) - uPrev(end))); % right boundary

    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    scalingW = ones(length(w),1);
    scalingW(1) = 0.5 * (1 + alf);

    kinEnergyW(n) = 1/2 * h * sum ((1/k * (w(2:end) - wPrev(2:end))).^2);
    connKinEnergyW(n) = 1/2 * h * scalingW(1) * (1/k * (w(1) - wPrev(1))).^2;
       
    potEnergyW(n) = c^2/(2 * h) * sum((w(2:end) - w(1:end-1)) .* (wPrev(2:end) - wPrev(1:end-1)));
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % left boundary

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

    connPotEnergy(n) = alf * c^2/(2 * h) * sum((interpolatedPoints(1) -  interpolatedPoints(2)) .* (interpolatedPointsPrev(1) - interpolatedPointsPrev(2)));
    totEnergy(n) = totEnergyU(n) + totEnergyW(n);% + connEnergy(n);
    totTotEnergy(n) = totEnergy(n) + connPotEnergy(n) + connKinEnergyU(n) + connKinEnergyW(n);
   
    %% save output
    outFree(n) = w(end - 5 * fs / 44100);

    %% draw stuff
    if n > drawStart && drawThings % && mod(n, drawSpeed) == 0
        
        gridMove = true;
        zoomed = true;
        addingPoint = true;
        if gridMove
            if addingPoint
                h = 1/31.5;
            else
                h = 1/30.5;
            end
            
            hLocsLeft = (0:(length(u))) * h * N;
            hLocsRight = (fliplr(1 - h * (0:(length(w))))) * N;
            
        else
            hLocsLeft = (0:(length(u)));
            hLocsRight = (fliplr(N - (0:(length(w)))));
            subplot(311)
        end
        hold off;
        uPlot = plot(hLocsLeft, [0;u], 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') ;

        hold on;
        wOffset = 0.00;
        if addingPoint
            wPlot = plot(hLocsRight(1:2), [w(1:2) + wOffset], 'Linewidth', 2, 'Color', 'r');
            scatter(hLocsRight(1), w(1), 80, 'r', 'Marker', 'o', 'Linewidth', 2);
            scatter(hLocsRight(2), w(2), 80, 'r', 'Marker', 'x', 'Linewidth', 2);
            plot(hLocsRight(2:3), [w(2:3) + wOffset], ':r', 'Linewidth', 2);
            scatter(hLocsRight(3), w(3), 50, 'r', 'Marker', 'o', 'Linewidth', 1);

        else   
            wPlot = plot(hLocsRight, [w + wOffset; 0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        end
        if gridMove
            xtickLocs = [0, hLocsLeft(ceil(length(u) / 2)), hLocsLeft(end), hLocsRight(1), hLocsRight(end-ceil(length(w) / 2)), N];
            xlabelSave = ["$u_0$", "$u_l$", "$u_M\ \ \;$", "$\ \ \;w_0$", "$w_l$", "$w_{M_w}$"];        
            if zoomed
                if ~addingPoint
                    scatter((hLocsRight(1) + hLocsRight(2)) * 0.5, 0, 40, 'b', 'Marker', 'o', 'Linewidth', 2)
                    scatter((hLocsLeft(end-1) + hLocsLeft(end)) * 0.5, 0, 300, 'r', 'Marker', '.', 'Linewidth', 2)
                    xtickLocs = [(hLocsLeft(end-1) + hLocsLeft(end)) * 0.5, hLocsLeft(end), hLocsRight(1), (hLocsRight(1) + hLocsRight(2)) * 0.5];
                    xlabelSave = ["$w_{-1}\ $", "$u_M$", "$w_0$", "$\ \ u_{N+1}$"];        
                else
                    scatter(hLocsLeft(end) + h * N, 0, 400, 'b', 'Marker', '.', 'Linewidth', 2)
                    plot([hLocsLeft(end), hLocsLeft(end) + h * N], [0, 0], 'b--', 'Linewidth', 2)
%                     xtickLocs = [hLocsLeft(end), hLocsLeft(end) + h * N, hLocsRight(1)];
%                     xlabelSave = ["$u_M$", "$I_3'\mathbf{v}^n$", "$w_0$"];        
                    xtickLocs = [hLocsLeft(end-1), hLocsLeft(end), hLocsRight(1), hLocsRight(2), hLocsRight(3)];
                    xlabelSave = ["$u_{M-1}$", "$u_M$", "$w_0$", "$w_1$", "$w_2$"];   
                    text(hLocsLeft(end) + h * N, 0.1, "$I_3'\mathbf{v}^n$", 'horizontalalignment', 'center', 'interpreter', 'latex', 'Fontsize', 16);
                end 
            end
        else
            xtickLocs = [0, floor(length(u) / 2), length(u), length(u) + floor(length(w) / 2), N];
            xlabelSave = ["$u_0$", "$u_l$", "$u_M, w_0$", "$w_l$", "$w_{M_w}$"];        
        end
        ylim([-0.6, 0.6])
        xlabel("$l$", 'interpreter', 'latex')
        grid on;
        if zoomed
            xlim([hLocsLeft(end-3), hLocsRight(4)])
        end
        ax = gca;
        ax.YTickLabel = ["$" + num2str(ax.YTick.') + repmat('\qquad',size(ax.YTickLabel,1),1) + "$"];
%         title("Sample = " + num2str(n) + "   N = " + num2str(floor(Ninit * 10) / 10))
        legend([uPlot, wPlot], ["$u$", "$w$"], 'Fontsize', 16, 'interpreter', 'latex')
        set(gca, 'Fontsize', 16, 'Linewidth', 2,...
            'xtick', xtickLocs, ...
            'xticklabel', xlabelSave, ...
            'TickLabelInterpreter', 'latex', ...
            'Position', [0.02 0.129186602870813 0.950000000000001 0.85645933014354]);
        set(gcf, 'Color', 'w');

        if n == 16
           disp("wait")
        end
        
%         subplot(312)
%         hold off;
%         plot(totEnergy(1:n) - totEnergy(1) + connKinEnergyU(1:n) + connKinEnergyW(1:n))
% 
%         hold on;
%         plot(connPotEnergy(1:n))
% 
%         subplot(313)
%         plot(totTotEnergy(1:n))

        drawnow;
%         if frame <= length(M) && filmFlag == true
%             M(frame) = getframe(gcf);
%             frame = frame + 1;
%         else
%             v = VideoWriter('dynamicGridCenter.mp4', 'MPEG-4');
%             v.FrameRate = 15;
%             open(v)
%             writeVideo(v, M);
%             close(v)
%             filmFlag = false;
%             break;
%         end
        
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
end

% subplot(2,1,1)
hold on
plot((1:lengthSound) / fs, outFree)
% 
% subplot(2,1,2)
% hold on;
% outfft = fft(outFree);
% plot([0:lengthSound-1]'*fs/lengthSound, 20 * log10(abs(outfft)), 'Linewidth', 2);
% xlim([0 0.4 * 44100])
% ylim([-60, 80])
% legend(["$f_s = 44100 \quad N =" + num2str(Ninit * 44100 / fs) + "$", ...
%     "$f_s =" + num2str(fs) + "\quad N =" + num2str(Ninit) + "$"], ...
%     'interpreter', 'latex')
% grid on
% xlabel("$f$ (in Hz)", 'interpreter', 'latex')
% ylabel("Magnitude (in dB)")
% set(gca, 'Fontsize', 16, 'Linewidth', 2)
% set(gcf, 'Color', 'w')
% 
% c/2