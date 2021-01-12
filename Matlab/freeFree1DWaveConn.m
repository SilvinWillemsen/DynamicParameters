clear all;
close all;
clc;

setting = false; % true is drawing, false is sound

drawSpeed = 1;
if drawSpeed == 1
    drawSpeedMod = 0;
else 
    drawSpeedMod = 1;
end
fs = 44100;             % Sample rate
k = 1/fs;               % Time step

lpConnection = true;
lpExponent = 30;

drawStart = 0;
drawThings = setting;
excite = true;

numFromBound = 1;

% if setting
    lengthSound = fs * 2;       % Length of the simulation
% else
%     lengthSound = fs * 1.5;
% end

startSample = 0;

changeC = true; % set to true for dynamic changes in wavespeed
sinusoidalChange = true;
freq = 15;

Ninit = 25.0 * fs / 44100;           % edit how many points you want
h = 1/Ninit;
Nend = 55.0 * fs / 44100;
cEnd = 1/(Nend*k);
cInit = h/k;            % calculate wave speed
c = cInit;

h = c*k;                % calculate h from wavespeed
N = floor(1/h);         % calculate points from h

lambdaSq = (c*k/h)^2    % should always be 1 as h is not recalculated

alf = Ninit - N;        % fractional remainder for the grid point
% alf = 0;
%% initialise states
if numFromBound == -1
    M = ceil(N/2);
    uNext = zeros(M, 1);
    u = zeros(M, 1);
else
    M = N-numFromBound;
    uNext = zeros(M, 1);
    u = zeros(M, 1);
end

excitationWidth = 0.2;
excitationLoc = 1/pi;
loc = excitationLoc * N;
width = max (4.0, excitationWidth * N);
raisedCosStart = floor (loc - width * 0.5);
raisedCosEnd = floor (loc + width * 0.5);

if excite
    u(raisedCosStart : raisedCosEnd) = 0.5 * (1 - cos (2.0 * pi * ((raisedCosStart:raisedCosEnd)-raisedCosStart) / width));
%     testPoint = floor(Ninit * 0.3 / 2);
%     u(floor(N/3.14)-testPoint:floor(N/3.14)+testPoint) = hann(testPoint * 2 + 1); % use hann window for excitation
%     u(end) = 1;
end
% u = rand(length(u), 1);
uPrev = u;


if numFromBound == -1
    wNext = zeros(N-M, 1);
    w = zeros(N-M, 1);
else
    wNext = zeros(numFromBound, 1);
    w = zeros(numFromBound, 1);
end
% w(floor(4*/5)-4:floor(4*N/5)+4) = hann(9); % use hann window for excitation

% w = rand(length(w), 1);
% w(1) = u(end);
wPrev = w;

% initialise laplacians
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

interpol = "sinc";
% custInterpol = interpol;
custInterpol = "sinc";
outFree = zeros(floor(lengthSound), 1);

if sinusoidalChange
    cVec = min(cInit, cEnd) + 0.5 * (abs(cEnd - cInit) * (1 + sign(cInit - cEnd) * cos (2 * pi * freq * (0:lengthSound-1) / fs)));
else
    cVec = linspace(cInit, cEnd, lengthSound);
end
%% recording
Mvid(500) = struct('cdata',[],'colormap',[]);
frame = 1;
filmFlag = true;
interpolatedPoints = [0; 0];

%% plotting
% figure('Position', [0, 0, 800, 200])
if setting
    figure('Position', [0, 0, 500, 200])
end

uSave = [];
wSave = [];
MuSave = [];
MwSave = [];

NVec = linspace(Ninit, Nend, lengthSound);
diffSave = zeros(lengthSound, 1);
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
  
    hLocsLeftTest = (0:(length(u))) * h;
    hLocsRightTest = (fliplr(1 - h * (0:(length(w)))));
    alfTest(n) = (hLocsRightTest(1) - hLocsLeftTest(end)) / h;
    
    if custInterpol == "cubic"
        customIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
                        2 * alf / (alf + 2); ...
                        2 / (alf + 2); ...
                        2 * alf / -((alf + 3) * (alf + 2))]';
        custSincWidth = 2;
    elseif custInterpol == "sinc"
        custSincWidth = 2;
        custBMax = pi;
        xPosCustIp = (-custSincWidth:custSincWidth-1) + [zeros(1, custSincWidth), ones(1, custSincWidth) * alf];
        customIp = sin(custBMax * xPosCustIp) ./ (xPosCustIp * custBMax);
        if alf == 0
            customIp(custSincWidth + 1) = 1;
        end
        
    end
    % add point if N^n > N^{n-1}
    if N ~= NPrev
        if N > NPrev
            if numFromBound == -1
                if mod(N,2) == 1
                    uNext = [uNext; customIp * [uNext(end-custSincWidth+1 : end); wNext(1:custSincWidth)]];
                    u = [u; customIp * [u(end-custSincWidth+1 : end); w(1:custSincWidth)]];
                    uPrev = [uPrev; customIp * [uPrev(end-custSincWidth+1 : end); wPrev(1:custSincWidth)]];
                else 
                    wNext = [fliplr(customIp) * [uNext(end-custSincWidth+1 : end); wNext(1:custSincWidth)]; wNext];
                    w = [fliplr(customIp) * [u(end-custSincWidth+1 : end); w(1:custSincWidth)]; w];
                    wPrev = [fliplr(customIp) * [uPrev(end-custSincWidth+1 : end); wPrev(1:custSincWidth)]; wPrev];
                end
            elseif numFromBound == 1
                uNext = [uNext; customIp(1:3) * [uNext(end-custSincWidth+1 : end); wNext(1)]];
                u = [u; customIp(1:3) * [u(end-custSincWidth+1 : end); w(1)]];
                uPrev = [uPrev; customIp(1:3) * [uPrev(end-custSincWidth+1 : end); wPrev(1)]];
            
            elseif numFromBound == 0
                disp("not defined");

            else
                uNext = [uNext; customIp * [uNext(end-custSincWidth+1 : end); wNext(1:custSincWidth)]];
                u = [u; customIp * [u(end-custSincWidth+1 : end); w(1:custSincWidth)]];
                uPrev = [uPrev; customIp * [uPrev(end-custSincWidth+1 : end); wPrev(1:custSincWidth)]];

            end
        end

        % remove point if N^n < N^{n-1}
        if N < NPrev
            diffSave(n) = w(1) - u(end);
            if numFromBound == -1 
                if mod(N,2) == 0
                    uNext = uNext(1:end-1);
                    u = u(1:end-1);
                    uPrev = uPrev(1:end-1);
                else 
                    wNext = wNext(2:end);
                    w = w(2:end);
                    wPrev = wPrev(2:end);
                end
            else
                uNext = uNext(1:end-1);
                u = u(1:end-1);
                uPrev = uPrev(1:end-1);
            end
            
        end
        eu = ones(length(u), 1);
        Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

        if numFromBound == -1
            M = ceil(N/2);
        else
            M = N-numFromBound;
        end
    end
    %% calculate interpolator
    if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
    elseif interpol == "sinc"
        includeUMp1AndWm1 = true;

        if alf < 1e-6
            alf = 1e-6;
%            includeUMp1AndWm1 = false;
        end
        alfSave(n) = alf;

%         sincWidth = floor(N / 2) - 1;
%         sincWidth = numFromBound + 1;
        sincWidth = 2;
        alphaBand = 1; % relative bandwidth range
        bmax = alphaBand*pi;

        if includeUMp1AndWm1
            xUMp1 = [-sincWidth:-1, -1:sincWidth-1]';
            xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth+1, 1)];
        else
            xUMp1 = (-sincWidth:sincWidth-1)';
            xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth, 1)];
        end
        iLen = length (xUMp1); % length of interpolation (N in stefans implementation)
        bU = (sin(bmax*xUMp1)./xUMp1);
        if sum(isnan(bU))
            idxIsNan = find(isnan(bU));
%                 bU(idxIsNan-2) = -bmax;
%                 bU(idxIsNan-1) = -bmax;
            bU(idxIsNan) = bmax;

        end 
        distU = xUMp1*ones(1,iLen)-ones(iLen,1)*xUMp1';    % distance matrix between points
        AU = sin(bmax*distU)./distU;
        AU(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
%             AU(isnan(AU)) = bmax;
        aU = AU\bU; %optimal coefficients
        
%         if alf < 1e-6
%             idxAu = find(round(aU) == 1);
%             aU(idxAu-2) = -1;
%             aU(idxAu-1) = 1;
%             aU(idxAu) = 1;
%         end

%% no need to recalculate for aW. We can just invert aU
%         if includeUMp1AndWm1
%             xWm1 = [-sincWidth+1:1, 1:sincWidth]';
%             xWm1 = xWm1 - [alf * ones(sincWidth+1, 1); zeros(sincWidth, 1)];
%         else
%             xWm1 = (-sincWidth+1:sincWidth)';
%             xWm1 = xWm1 - [alf * ones(sincWidth, 1); zeros(sincWidth, 1)];
%         end
%         bW = (sin(bmax*xWm1)./xWm1);
%         if sum(isnan(bW))
%             bW(isnan(bW)) = bmax;
%         end 
%         distW = xWm1*ones(1,iLen)-ones(iLen,1)*xWm1';    % distance matrix between points
%         AW = sin(bmax*distW)./distW;
%         AW(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
% %             AW(isnan(AW)) = bmax;
%         aW = AW\bW; %optimal coefficients
%         
        aW = flipud(aU);

        sincRange = (M-sincWidth : M+sincWidth);
    else
        ip = [0, (1-alf), alf, 0];
    end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    % calculate interpolated points
    interpolatedPointsPrev = interpolatedPoints;

    if lpConnection
        diffAtConn = w(1) - u(end);
        diffAtConnPrev = wPrev(1) - uPrev(end);
        u(end) = u(end) + (1-alf)^lpExponent * diffAtConn * 0.5;
        w(1) = w(1) - (1-alf)^lpExponent * diffAtConn * 0.5;
    end
    
    if interpol == "cubic"
        if numFromBound == 1
            interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [(ip(3) - ip(1)) * w(1);
                                            ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
        else
            interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                            ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
        end
    elseif interpol == "sinc"
        curState = [u; w];
        if M + sincWidth == N+1 % use boundary condition
            if includeUMp1AndWm1
                interpolatedPoints = [(aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)]) * curState(M-sincWidth+1:end);...
                                       aW(1:end-1)' * curState(M-sincWidth : end)];
            else
                interpolatedPoints = [(aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)]) * curState(M-sincWidth+2:end);...
                       aW(1:end-1)' * curState(M-sincWidth : end-1)];

            end
%             BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW(1:end-1)';
        elseif M + sincWidth == N % use boundary condition
            if includeUMp1AndWm1
                interpolatedPoints = [aU(1:end-1)' * curState(M-sincWidth+1:end);...
                                       aW' * curState(M-sincWidth : end)];
            else
                interpolatedPoints = [aU(1:end-1)' * curState(M-sincWidth+2:end);...
                       aW' * curState(M-sincWidth : end-1)];

            end
%             interpolatedPoints = [aU(1:end-1)' * curState(M-sincWidth+1:end); aW' * curState(M-sincWidth : end)];
        else
            % still need to exclude UM and W0 if ~includeUMp1AndWm1
            if includeUMp1AndWm1
                sincRange = (M-sincWidth : M+sincWidth);
            else
                sincRange =[M-sincWidth:M-1, M+1:M+sincWidth];
            end
            interpolatedPoints = [aU' * curState(sincRange+1); aW' * curState(sincRange)];
        end
    end
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    uNext(end) = uNext(end) + lambdaSq * interpolatedPoints(1);
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + lambdaSq * interpolatedPoints(2);

    
%     wIp = wNext(1);
%     uIp = uNext(end);
%     diffAtConn = wIp - uIp;
%     uNext(end) = uNext(end) + (1-alf)^10 * diffAtConn * 0.5;
%     wNext(1) = wNext(1) - (1-alf)^10 * diffAtConn * 0.5;
%     
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
%     kinEnergyU(n) = 1/2 * h * sum ((1/k * (u(1:end-1) - uPrev(1:end-1))).^2);
%     connKinEnergyU(n) = 1/2 * h * scalingU(end) * (1/k * (u(end) - uPrev(end))).^2;
%     
%     potEnergyU(n) = c^2/(2 * h) * sum((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)));
%     potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((u(1) - 0) .* (uPrev(1) - 0)); % left boundary
% %     potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * sum((w(2) - u(end)) .* (wPrev(2) - uPrev(end))); % right boundary
% 
%     totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
%     
%     scalingW = ones(length(w),1);
%     scalingW(1) = 0.5 * (1 + alf);
% 
%     kinEnergyW(n) = 1/2 * h * sum ((1/k * (w(2:end) - wPrev(2:end))).^2);
%     connKinEnergyW(n) = 1/2 * h * scalingW(1) * (1/k * (w(1) - wPrev(1))).^2;
%        
%     potEnergyW(n) = c^2/(2 * h) * sum((w(2:end) - w(1:end-1)) .* (wPrev(2:end) - wPrev(1:end-1)));
%     potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % left boundary
% 
%     totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);
% 
%     connPotEnergy(n) = alf * c^2/(2 * h) * sum((interpolatedPoints(1) -  interpolatedPoints(2)) .* (interpolatedPointsPrev(1) - interpolatedPointsPrev(2)));
%     totEnergy(n) = totEnergyU(n) + totEnergyW(n);% + connEnergy(n);
%     totTotEnergy(n) = totEnergy(n) + connPotEnergy(n) + connKinEnergyU(n) + connKinEnergyW(n);
%    
    %% save output
    outFree(n) = u (fs / 44100);

    %% draw stuff
    if n > drawStart && drawThings && mod(n, drawSpeed) == drawSpeedMod
        
        gridMove = true;
        zoomed = false;       
        interpolExplanation = false;

        addingPoint = false;
        if gridMove
%             if addingPoint
%                 h = 1/31.5;
%             else
%                 h = 1/30.5;
%             end
            
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
            if numFromBound ~= 1
                wPlot = plot(hLocsRight(1:2), [w(1:2) + wOffset], 'Linewidth', 2, 'Color', 'r');
                scatter(hLocsRight(1), w(1), 80, 'r', 'Marker', 'o', 'Linewidth', 2);
                scatter(hLocsRight(2), w(2), 80, 'r', 'Marker', 'x', 'Linewidth', 2);
                plot(hLocsRight(2:3), [w(2:3) + wOffset], ':r', 'Linewidth', 2);
                scatter(hLocsRight(3), w(3), 50, 'r', 'Marker', 'o', 'Linewidth', 1);
            else
                wPlot = plot([hLocsRight(1); hLocsRight(1) + 1], [w(1); 0] + wOffset, 'Linewidth', 2, 'Color', 'r');
                scatter(hLocsRight(1), w(1), 80, 'r', 'Marker', 'o', 'Linewidth', 2);
                scatter(hLocsRight(1) + 1, 0, 80, 'r', 'Marker', 'x', 'Linewidth', 2);
                plot([hLocsRight(1) + 1; hLocsRight(1) + 2], [0; -w(1)] + wOffset, ':r', 'Linewidth', 2);
                scatter(hLocsRight(1) + 2, -w(1), 50, 'r', 'Marker', 'o', 'Linewidth', 1);
            end

        else   
            wPlot = plot(hLocsRight, [w + wOffset; 0], 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        end
        if gridMove
            eps = 1e-6;
            xtickLocs = [0, hLocsLeft(ceil(length(u) / 2)), hLocsLeft(end), hLocsRight(1) + eps, hLocsRight(end-ceil(length(w) / 2)), N];
            xlabelSave = ["$u_0$", "$u_l$", "$u_M\ \ \;$", "$\ \ \;w_0$", "$w_l$", "$w_{M_w}$"];        
            if zoomed
                if ~addingPoint
                    scatter(hLocsLeft(end) + 1, interpolatedPoints(1), 40, 'b', 'Marker', 'o', 'Linewidth', 2)
                    scatter(hLocsRight(1) - 1, interpolatedPoints(2), 300, 'r', 'Marker', '.', 'Linewidth', 2)
                    if hLocsLeft(end) == hLocsRight(1)
                        xtickLocs = [hLocsRight(1) - 1, hLocsLeft(end), hLocsLeft(end) + 1];
                        xlabelSave = ["$w_{-1}\ $", "$u_M, w_0$", "$\ \ u_{M+1}$"]; 
                    else 
                        xtickLocs = [hLocsRight(1) - 1, hLocsLeft(end), hLocsRight(1), hLocsLeft(end) + 1];
                        xlabelSave = ["$w_{-1}\ $", "$u_M$", "$w_0$", "$\ \ u_{M+1}$"]; 
                    end
                else
                    scatter(hLocsLeft(end) + 1, 0, 40, 'b', 'Marker', 'o', 'Linewidth', 2)
                    scatter(hLocsRight(1) - 1, 0, 400, 'r', 'Marker', '.', 'Linewidth', 2)

%                     plot([hLocsLeft(end), hLocsLeft(end) + h * N], [0, 0], 'b--', 'Linewidth', 2)
%                     xtickLocs = [hLocsLeft(end), hLocsLeft(end) + h * N, hLocsRight(1)];
%                     xlabelSave = ["$u_M$", "$I_3'\mathbf{v}^n$", "$w_0$"];        
                     if hLocsLeft(end) == hLocsRight(1)
                        xtickLocs = [hLocsRight(1) - 1, hLocsLeft(end), hLocsLeft(end) + 1, hLocsRight(1) + 2];
                        xlabelSave = ["$w_{-1}, u_{M-1}$", "$u_M, w_0$", "$\ u_{M+1}, w_1$","$w_2$"]; 
                    else 
                        xtickLocs = [hLocsLeft(end-1), hLocsRight(1) - 1, hLocsLeft(end), hLocsRight(1), hLocsLeft(end) + 1, hLocsRight(1) + 1, hLocsRight(1) + 2];
                        xlabelSave = ["$u_{M-1}$", "$w_{-1}\ $", "$u_M$", "$w_0$", "$\ \ u_{M+1}$", "$w_1$", "$w_2$"]; 
                     end
%                      text(hLocsLeft(end) + h * N, 0.1, "$\mathcal{I}\mathbf{v}^n$", 'horizontalalignment', 'center', 'interpreter', 'latex', 'Fontsize', 16);
                    grid on;
                    yticks([0])
                    text((hLocsRight(1) + hLocsLeft(end)) * 0.5, 0.1, "$\alpha = "+ num2str(round(alf * 100) / 100) + "$", ...
                        'interpreter', 'latex', 'Fontsize', 16, ...
                        'horizontalAlignment', 'center', ...
                        'color', [0.5, 0.5, 0.5]);
                end 
            end
        else
            xtickLocs = [0, floor(length(u) / 2), length(u), length(u) + floor(length(w) / 2), N];
            xlabelSave = ["$u_0$", "$u_l$", "$u_M, w_0$", "$w_l$", "$w_{M_w}$"];        
        end
%         ylim([-0.6, 0.6])
%         ylim([-5,5])
        xlabel("$l$", 'interpreter', 'latex')
        grid on;
        if zoomed
%             xlim([hLocsLeft(end-3), hLocsRight(4)])
            if numFromBound == -1
%                 xlim([hLocsLeft(12), hLocsRight(5)])
                xlim([0.45 * N, 0.55 * N])
            else
                xlim([N-numFromBound - N/4, N-numFromBound + N/8])
            end
        else
            xlim([0, N])
        end
        if interpolExplanation
%             xlim([hLocsLeft(14), hLocsLeft(14) + 6])

            alfOffset = -1.2;
            text((hLocsLeft(end) + hLocsRight(1)) * 0.5, alfOffset, "$\alpha =" + num2str(alf, 2) +"$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 14, 'color', [0.5, 0.5, 0.5]);
            plot([hLocsLeft(end), hLocsLeft(end)], [alfOffset*0.8, u(end)], '--', 'color', [0.5, 0.5, 0.5])
            plot([hLocsRight(1), hLocsRight(1)], [alfOffset*0.8, u(end)], '--', 'color', [0.5, 0.5, 0.5])
            plot([hLocsLeft(end), hLocsRight(1)], [alfOffset*0.8, alfOffset*0.8], 'color', [0.5, 0.5, 0.5])
            grid off
            axis off
            offset = 0.3;
            text(hLocsLeft(end-1), u(end-1) - offset, "$u_{M-1}^n$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 16, 'color', 'b');
            text(hLocsLeft(end), u(end) - offset, "$u_{M}^n$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 16, 'color', 'b');
            text(hLocsLeft(end)+1, interpolatedPoints(1) - offset, "$u_{M+1}^n$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 16, 'color', 'b');
            
            text(hLocsRight(1) - 1, interpolatedPoints(2) + offset, "$w_{-1}^n$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 16, 'color', 'r');
            text(hLocsRight(1), w(1) + offset, "$w_0^n$", ... 
                'interpreter', 'latex', 'horizontalAlignment', 'center', ...
                'Fontsize', 16, 'color', 'r');
%             text(hLocsRight(2), w(2) + offset, "$w_1^n$", ... 
%                 'interpreter', 'latex', 'horizontalAlignment', 'center', ...
%                 'Fontsize', 16, 'color', 'r');
            set(gcf, 'color', 'w')
            ylim([-1.5, 1.5])
        end
        ax = gca;
        ax.YTickLabel = ["$" + num2str(ax.YTick.') + repmat('\qquad',size(ax.YTickLabel,1),1) + "$"];
%         title("Sample = " + num2str(n) + "   N = " + num2str(floor(Ninit * 10) / 10))
        legend([uPlot, wPlot], ["$u$", "$w$"], 'Fontsize', 16, 'interpreter', 'latex')
        set(gca, 'Fontsize', 16, 'Linewidth', 2,...
            'TickLabelInterpreter', 'latex', ...
            'xticklabel', xlabelSave, ... % 'XTick', xtickLocs, ...
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
if ~setting
    figure
    spectrogram(outFree,512,64,512, fs)
    view(90, -90)
    set(gcf, 'color', 'w')
end
% subplot(2,1,1)
hold on
% plot((1:lengthSound) / fs, outFree)
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