clear all;
close all;
clc;

setting = true; % true is drawing, false is sound
plotEnergy = true;
plotROCEnergy = false;
drawSpeed = 1;
if drawSpeed == 1
    drawSpeedMod = 0;
else 
    drawSpeedMod = 1;
end
fs = 44100;             % Sample rate
k = 1/fs;               % Time step

dispCorr = false;
corrAlf = false;
lpConnection = false;
lpExponent = 30;

drawStart = 0;
drawThings = setting;
excite = true;
exciteConnection = ~excite;

numFromBound = -1;

%{
Decides whether to use the full range or not 
    0: fullSinc is false 
    1: include all moving points
    2: include boundaries as well
    3: include virtual grid points
%}
fullSinc = 3; 
fSCentered = false; % "true" only works if numFromBound == -1

% if setting
    lengthSound = fs;       % Length of the simulation
% else
%     lengthSound = fs * 1.5;
% end

startSample = 0;

changeC = true; % set to true for dynamic changes in wavespeed
sinusoidalChange = false;
freq = 10;

Ninit = 160.0 * fs / 44100;           % edit how many points you want
h = 1/Ninit;
Nend = 160.0 * fs / 44100;
cEnd = 1/(Nend*k);
cInit = h/k;            % calculate wave speed
c = cInit;

alf = Ninit - floor(Ninit);        % fractional remainder for the grid point
% sigma1 = 0.0005;
% sig = sigma1 / ((alf + 1) * 0.5);
hScalar = 1;
% h = hScalar * sqrt((c*k)^2 + 2 * sig * k);                % calculate h from wavespeed
h = c * k; 
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

excitationWidth = 0.1;
excitationLoc = 4/5;
loc = excitationLoc * N;
width = floor(max (2.0, excitationWidth * N));
raisedCosStart = floor (loc - width * 0.5);
raisedCosEnd = raisedCosStart + width;
exciteForFigure = false;

if excite && loc < M
        u(raisedCosStart : raisedCosEnd) = 0.5 * (1 - cos (2.0 * pi * ((raisedCosStart:raisedCosEnd)-raisedCosStart) / width));
        %     testPoint = floor(Ninit * 0.3 / 2);
%     u(floor(N/3.14)-testPoint:floor(N/3.14)+testPoint) = hann(testPoint * 2 + 1); % use hann window for excitation
%     u(end-1) = 1;
end

if exciteConnection
    u(end) = 0.1;
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

if excite && loc > M
    raisedCosStart = raisedCosStart - M;
    raisedCosEnd = raisedCosEnd - M;
    w(raisedCosStart : raisedCosEnd) = 0.5 * (1 - cos (2.0 * pi * ((raisedCosStart:raisedCosEnd)-raisedCosStart) / width));
end


if exciteForFigure
    if Ninit ~= 30
        disp("not correct number of points")
    else
        u = zeros(size(u));
        w = zeros(size(u));
        w(end-8:end) = hann(9);
    end
end

% w(floor(4*/5)-4:floor(4*N/5)+4) = hann(9); % use hann window for excitation

% w = rand(length(w), 1);
% w(1) = u(end);
if exciteConnection
    w(1) = -0.1;
end
wPrev = w;

% u = uSave;
% uPrev = uPrevSave;
% w = wSave;
% wPrev = wPrevSave;
% initialise laplacians
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

interpol = "quadratic";
% custInterpol = interpol;
custInterpol = "cubic";
outFree = zeros(floor(lengthSound), 1);

if sinusoidalChange
    cVec = min(cInit, cEnd) + 0.5 * (abs(cEnd - cInit) * (1 + sign(cInit - cEnd) * cos (2 * pi * freq * (0:lengthSound-1) / fs)));
else
%     cVec = [linspace(cInit, cEnd, floor(lengthSound/2)), cEnd * ones(1, ceil(lengthSound/2))];
    cVec = linspace(cInit, cEnd, floor(lengthSound));

end
%% recording
Mvid(200) = struct('cdata',[],'colormap',[]);
frame = 1;
filmFlag = false;
interpolatedPoints = [0; 0];

%% plotting
% figure('Position', [0, 0, 800, 200])
if setting
    if plotEnergy
        figure('Position', [0, 0, 700, 500])
    else
        figure('Position', [0, 0, 500, 130])
    end
end

uSave = [];
wSave = [];
MuSave = [];
MwSave = [];

NVec = linspace(Ninit, Nend, lengthSound);
diffSave = zeros(lengthSound, 1);
diffSaveIdx = 0;

boundaryEnergyUH = zeros(lengthSound, 1);
boundaryEnergyWH = zeros(lengthSound, 1);

for n = 1:lengthSound  

    % change wave speed
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
%     sig = sigma1 / ((alf + 1) / 2);
    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
%     h = hScalar * sqrt((c*k)^2 + 2 * sig * k);
%     sigSave(n) = sig;

    h = c * k;
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
            diffSaveIdx = diffSaveIdx + 1;
            diffSave(diffSaveIdx) = w(1) - u(end);
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
    if interpol == "quadratic"
        ip = [(alf - 1) / (alf + 1), 1, -(alf - 1) / (alf + 1)];
    elseif interpol == "cubic"
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

        if numFromBound == -1
            sincWidth = floor(N / 2);
        else
            sincWidth = numFromBound + 1;
        end
        sincWidth = 2;
        alphaBand = 1; % relative bandwidth range
        bmax = alphaBand*pi;

        if fullSinc == 0
            if includeUMp1AndWm1
                xUMp1 = [-sincWidth:-1, -1:sincWidth-1]';
                xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth+1, 1)];
            else
                xUMp1 = (-sincWidth:sincWidth-1)';
                xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth, 1)];
            end
        elseif fullSinc == 1
            xUMp1 = (1:N)' - M - 1;
            xUMp1 = xUMp1 + [zeros(M, 1); alf * ones(N-M, 1) - 1];
        elseif fullSinc == 2
            xUMp1 = (0:N+1)' - M - 1;
            xUMp1 = xUMp1 + [zeros(M+1, 1); alf * ones(N-M+1, 1) - 1];
        elseif fullSinc == 3
            xUMp1 = (-1:N+2)' - M - 1;
            xUMp1 = xUMp1 + [zeros(M+2, 1); alf * ones(N-M+2, 1) - 1];
        end
        if numFromBound == -1 && fSCentered && fullSinc ~= 0
            xUMp1 = xUMp1(2:end);
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
        if includeUMp1AndWm1
            xWm1 = [-sincWidth+1:1, 1:sincWidth]';
            xWm1 = xWm1 - [alf * ones(sincWidth+1, 1); zeros(sincWidth, 1)];
        else
            xWm1 = (-sincWidth+1:sincWidth)';
            xWm1 = xWm1 - [alf * ones(sincWidth, 1); zeros(sincWidth, 1)];
        end
        if fullSinc == 0
%             aW = flipud(aU);
            sincRange = (M-sincWidth : M+sincWidth);
        elseif fullSinc == 1
            xWm1 = (1:N)' - M;
            xWm1 = xWm1 - [alf * ones(M, 1) - 1; zeros(N-M, 1)];
        elseif fullSinc == 2
            xWm1 = (0:N+1)' - M;
            xWm1 = xWm1 - [alf * ones(M+1, 1) - 1; zeros(N-M+1, 1)];

        elseif fullSinc == 3
            xWm1 = (-1:N+2)' - M;
            xWm1 = xWm1 - [alf * ones(M+2, 1) - 1; zeros(N-M+2, 1)];
        end
        if numFromBound == -1 && fSCentered && fullSinc ~= 0
            xWm1 = xWm1(1:end-1);
        end
%         if fullSinc ~= 0
        bW = (sin(bmax*xWm1)./xWm1);
        if sum(isnan(bW))
            bW(isnan(bW)) = bmax;
        end
        distW = xWm1*ones(1,iLen)-ones(iLen,1)*xWm1';    % distance matrix between points
        AW = sin(bmax*distW)./distW;
        AW(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
        AW(isnan(AW)) = bmax;
        aW = AW\bW; %optimal coefficients
%         end

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
        lpVec = 0.5 * diffAtConn * [-(1-alf)^lpExponent, (1-alf)^lpExponent];
%         lpVecPrev = 0.5 * diffAtConnPrev * [-(1-alf)^lpExponent, (1-alf)^lpExponent];

%         lpVec = 0.5 * diffAtConn * [-cos((alf) * pi/2)^lpExponent, cos((alf) * pi/2)^lpExponent];

%         u(end) = u(end) + (1-alf)^lpExponent * diffAtConn * 0.5;
%         w(1) = w(1) - (1-alf)^lpExponent * diffAtConn * 0.5;
        
        u(end) = u(end) + lpVec(2);
        w(1) = w(1) + lpVec(1);
%         
%         uPrev(end) = uPrev(end) + lpVecPrev(2);
%         wPrev(1) = wPrev(1) + lpVecPrev(1);
%         uDampTerm = sig * k / (h * h) * (w(1) - 2 * u(end) + u(end-1) ... 
%             - wPrev(1) + 2 * uPrev(end) - uPrev(end-1));
%         wDampTerm = sig * k / (h * h) * (w(2) - 2 * w(1) + u(end) ... 
%             - wPrev(2) + 2 * wPrev(1) - uPrev(end-1));
    end
    
    if interpol == "quadratic"
        if numFromBound == 1            
            interpolatedPoints = [ip(1) * u(end) + ip(2) * w(1) + 0; ...
                                  ip(1) * w(1) + ip(2) * u(end) + ip(3) * u(end-1)];
        else
            interpolatedPoints = [ip(1) * u(end) + ip(2) * w(1) + ip(3) * w(2); ...
                                  ip(1) * w(1) + ip(2) * u(end) + ip(3) * u(end-1)];
        
        end
    elseif interpol == "cubic"
        if numFromBound == 1
            interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [(ip(3) - ip(1)) * w(1);
                                            ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
        else
            interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                            ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
        end
    elseif interpol == "sinc"
        if numFromBound == -1 && fSCentered
            inputRange1 = 2:N;
            inputRange2 = 1:N-1;
        else
            inputRange1 = 1:N;
            inputRange2 = 1:N;
        end
        curState = [u; w];
        if fullSinc == 1
            interpolatedPoints = [aU' * curState(inputRange1); aW' * curState(inputRange2)];
        elseif fullSinc == 2
            interpolatedPoints = [aU(2:end-1)' * curState(inputRange1); aW(2:end-1)' * curState(inputRange2)];
        elseif fullSinc == 3
            interpolatedPoints = [(aU(3:end-2)' - [aU(1), zeros(1, length(aU) - 6), aU(end)]) * curState(inputRange1);...
                (aW(3:end-2)' - [aW(1), zeros(1, length(aW) - 6), aW(end)]) * curState(inputRange2)];
        elseif M + sincWidth == N+1 % use boundary condition
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
%     if lpConnection
%         uNext(end) = uNext(end) + uDampTerm;
%     end
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + lambdaSq * interpolatedPoints(2);
    
    if dispCorr
        epsilon = 0;

        etaDiv = 1;
        
        eta = (w(1) - u(end)) * etaDiv;
        etaPrev = (wPrev(1) - uPrev(end)) * etaDiv;
               
        if corrAlf && alf < 0.001
            sig0 = k;
            disp("alf < 0.001")
%             if n ~= 1
%                 drawThings = true;
%             end
        else
            sig0 = 1;
        end
        rForce = (1 - sig0 / k) / (1 + sig0 / k);
        oOP = (h * (1 + sig0 / k) * (1-alf)) / (2 * h * (alf + epsilon) + 2 * etaDiv * k^2 * (1 + sig0 / k) * (1-alf));
        
        F = ((wNext(1) - uNext(end)) * etaDiv + rForce * etaPrev) * oOP;
        
%         oOP = ((1-alf) * h * k + 2 * sig0 * h * (alf + epsilon)) / (2 * (alf + epsilon) * k * h + (1-alf) * k^3 + 2 * sig0 * k^2 * (alf + epsilon));
%         rForce = ((1-alf) * k - 2 * sig0 * (alf + epsilon)) / ((1-alf) * k + 2 * sig0 * (alf + epsilon));
%         
%         F = ((wNext(1) - uNext(end)) * etaDiv + rForce * etaPrev) * oOP;

%         F = ((wNext(1) - uNext(end)) * etaDiv + etaPrev) * h / (2 * alf * h + 2 * epsilon * h + k^2);
%         F = ((wNext(1) - uNext(end)) * etaDiv + etaPrev) * (1-alf) *  h / ((2 * alf + epsilon) * h + k^2 * (1-alf));
               
%         F = (wNext(1) - uNext(end) + 4 * eta + 2 * etaPrev) * h / (8 * alf * h^2 + 8 * epsilon * h + 2 * k^2);

        uNext(end) = uNext(end) + k^2/h * F;
        wNext(1) = wNext(1) - k^2/h * F;

    end
    
%     if lpConnection
%         wNext(1) = wNext(1) + wDampTerm;
%     end
    
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
    epsilonUr = 1 + alf;
    epsilonWr = 1 + alf;

    
    scalingU = ones(length(u),1);
    scalingU(end) = 0.5 * (1 + alf);
    kinEnergyU(n) = 1/2 * h * sum (scalingU .* (1/k * (u - uPrev)).^2);
%     kinEnergyU(n) = 1/2 * h * sum ((1/k * (u(1:end-1) - uPrev(1:end-1))).^2);
%     connKinEnergyU(n) = 1/2 * h * scalingU(end) * (1/k * (u(end) - uPrev(end))).^2;
    
    potEnergyU(n) = c^2/(2 * h) * sum((u(2:end) - u(1:end-1)) .* (uPrev(2:end) - uPrev(1:end-1)));
    potEnergyU(n) = potEnergyU(n) + c^2/(2 * h) * (u(1) - 0) * (uPrev(1) - 0); % left boundary

    idx = n - (1  * (n>1));
    boundaryEnergyU(n) = c^2 * epsilonUr / (16 * h * k) * (uNext(end) - uPrev(end)) * (interpolatedPoints(1) - 2 * u(end) + u(end-1));
    boundaryEnergyUH(n) = k * boundaryEnergyU(n) + boundaryEnergyUH(idx);
    
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    scalingW = ones(length(w),1);
    scalingW(1) = 0.5 * (1 + alf);
    
    kinEnergyW(n) = 1/2 * h * sum ((1/k * scalingW .* (w - wPrev)).^2);

%     kinEnergyW(n) = 1/2 * h * sum ((1/k * (w(2:end) - wPrev(2:end))).^2);
%     connKinEnergyW(n) = 1/2 * h * scalingW(1) * (1/k * (w(1) - wPrev(1))).^2;
       
    potEnergyW(n) = c^2/(2 * h) * sum((w(2:end) - w(1:end-1)) .* (wPrev(2:end) - wPrev(1:end-1)));
    potEnergyW(n) = potEnergyW(n) + c^2/(2 * h) * sum((0 - w(end)) .* (0 - wPrev(end))); % right boundary
    
    boundaryEnergyW(n) = c^2 * epsilonWr / (16 * h * k) * (wNext(1) - wPrev(1)) * (interpolatedPoints(2) - 2 * w(1) + w(2));
    boundaryEnergyWH(n) = k * boundaryEnergyW(n) + boundaryEnergyWH(idx);

    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

%     connPotEnergyU(n) = c^2/(2 * h) * (interpolatedPoints(1) -  u(end)) * (interpolatedPointsPrev(1) - uPrev(end));
%     connPotEnergyW(n) = c^2/(2 * h) * (w(1) -  interpolatedPoints(2)) * (wPrev(1) - interpolatedPointsPrev(2));

%     connPotEnergyU(n) = c^2/(2 * h) * (w(2) -  u(end)) * (wPrev(2) - uPrev(end));
%     connPotEnergyW(n) = c^2/(2 * h) * (w(1) -  u(end-1)) * (wPrev(1) - uPrev(end-1));


    totEnergy(n) = totEnergyU(n) + totEnergyW(n);% + connEnergy(n);
    totTotEnergy(n) = totEnergy(n) + boundaryEnergyUH(idx) + boundaryEnergyWH(idx);% + connKinEnergyU(n) + connKinEnergyW(n);% + connPotEnergyU(n) + connPotEnergyU(n);
%   

    %% ROC energies
    Mw = length(w);

    scalingROCU = ones(M,1);
    scalingROCU(end) = epsilonUr / 2;

    scalingROCW = ones(length(w),1);
    scalingROCW(1) = epsilonWr / 2;
    
    potEnergyRangeU = 2:M-1;
    potEnergyRangeW = 2:Mw-1;

    rOCKinEnergyU(n) = h / (2 * k^3) * sum(scalingROCU .* (uNext - uPrev) .* (uNext - 2 * u + uPrev));
    rOCPotEnergyU(n) = c^2 / (2 * h * k) * sum((uNext(potEnergyRangeU) - uPrev(potEnergyRangeU)) ...
        .* (u(potEnergyRangeU+1) - 2 * u(potEnergyRangeU) + u(potEnergyRangeU-1))) ...
        + c^2 * epsilonUr / (4 * h * k) * (uNext(M) - uPrev(M)) * (interpolatedPoints(1) - 2 * u(M) + u(M-1));
    rOCconnectionU(n) = c^2 * epsilonUr / (4 * h * k) * (uNext(M) - uPrev(M)) * (interpolatedPoints(1) - 2 * u(M) + u(M-1));
    
    rOCPotEnergyU(n) = rOCPotEnergyU(n) + c^2 / (2 * h * k) * (uNext(1) - uPrev(1)) ...
        * (u(2) - 2 * u(1) + 0); % left boundary
    
    rOCKinEnergyW(n) = h / (2 * k^3) * sum(scalingROCW .* (wNext - wPrev) .* (wNext - 2 * w + wPrev));
    rOCPotEnergyW(n) = c^2 / (2 * h * k) * sum((wNext(potEnergyRangeW) - wPrev(potEnergyRangeW)) ...
        .* (w(potEnergyRangeW+1) - 2 * w(potEnergyRangeW) + w(potEnergyRangeW-1))) ...
        + c^2 * epsilonWr / (4 * h * k) * (wNext(1) - wPrev(1)) * (interpolatedPoints(2) - 2 * w(1) + w(2));
    rOCconnectionW(n) = c^2 * epsilonWr / (4 * h * k) * (wNext(1) - wPrev(1)) * (interpolatedPoints(2) - 2 * w(1) + w(2));
    rOCPotEnergyW(n) = rOCPotEnergyW(n) + c^2 / (2 * h * k) * (wNext(end) - wPrev(end)) ...
        * (0 - 2 * w(end) + w(end-1)); % right boundary

    totROCenergy(n) = rOCKinEnergyU(n) - rOCPotEnergyU(n) + rOCKinEnergyW(n) - rOCPotEnergyW(n);% + rOCconnectionU(n) - rOCconnectionW(n);
    
    

    %% save output
    outFree(n) = u (fs / 44100);

    %% draw stuff
    if n > drawStart && drawThings && mod(n, drawSpeed) == drawSpeedMod
        if ~plotEnergy
            gridMove = true;
            zoomed = false;  
            lpExplanation = false;
            interpolExplanation = false;
            includeAlpha  = false;
            addingPoint = false;
            
            if lpExplanation
                u(end) = -0.1;
                w(1) = 0.1;
            end
            if gridMove
%                 if addingPoint
%                     h = 1/31.5;
%                 else
%                     h = 1/30.5;
%                 end

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
                    wPlot = plot(hLocsRight(1:end-1), [w + wOffset], 'Linewidth', 2, 'Color', 'r');
                    scatter(hLocsRight(1:end-1), w, 80, 'r', 'Marker', 'o', 'Linewidth', 2);
    %                 scatter(hLocsRight(2), w(2), 80, 'r', 'Marker', 'x', 'Linewidth', 2);
    %                 plot(hLocsRight(2:3), [w(2:3) + wOffset], ':r', 'Linewidth', 2);
    %                 scatter(hLocsRight(3), w(3), 50, 'r', 'Marker', 'o', 'Linewidth', 1);
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
                        if ~lpExplanation
                            scatter(hLocsLeft(end) + 1, interpolatedPoints(1), 40, 'b', 'Marker', 'o', 'Linewidth', 2)
                            scatter(hLocsRight(1) - 1, interpolatedPoints(2), 300, 'r', 'Marker', '.', 'Linewidth', 2)
                        end
                        if hLocsLeft(end) == hLocsRight(1)
                            if lpExplanation
                                xtickLocs = [hLocsLeft(end-1), hLocsLeft(end), hLocsRight(2)];
                                xlabelSave = ["$u_{M-1}\ $", "$u_M, w_0$", "$\ \ w_{1}$"]; 
                            else
                                xtickLocs = [hLocsRight(1) - 1, hLocsLeft(end), hLocsLeft(end) + 1];
                                xlabelSave = ["$w_{-1}\ $", "$u_M, w_0$", "$\ \ u_{M+1}$"]; 
                            end
                        else 
                            if lpExplanation
                                xtickLocs = [hLocsLeft(end-1), hLocsLeft(end), hLocsRight(1), hLocsRight(2)];
                                xlabelSave = ["$u_{M-1}\ $", "$u_M$", "$w_0$", "$\ \ w_{1}$"]; 
                            else
                                xtickLocs = [hLocsLeft(end-1), hLocsRight(1) - 1, hLocsLeft(end), hLocsRight(1), hLocsLeft(end) + 1, hLocsRight(2)];
                                xlabelSave = ["$u_{M-1}\ \ \ $", "$w_{-1}$", "$u_M$", "$w_0$", "$u_{M+1}$", "\ $w_1$"];
                                if includeAlpha                                    
                                    lHeight = -0.15;
                                    lWidth = 0.025;
                                    text((hLocsLeft(end) + hLocsRight(1)) / 2, lHeight + sign(lHeight) * 0.15, "$\alpha = 0.5$", 'horizontalalignment', 'center',...
                                        'interpreter', 'latex', 'Fontsize', 16, ...
                                        'color', [0.4, 0.4, 0.4]);

                                    plot ([hLocsLeft(end), hLocsRight(1)], [lHeight, lHeight], 'Linewidth', 1,  'color', [0.4, 0.4, 0.4]);
                                    plot ([hLocsLeft(end), hLocsLeft(end)], [lHeight + lWidth, lHeight - lWidth], 'Linewidth', 1,  'color', [0.4, 0.4, 0.4]);
                                    plot ([hLocsRight(1), hLocsRight(1)], [lHeight + lWidth, lHeight - lWidth], 'Linewidth', 1,  'color', [0.4, 0.4, 0.4]);
%                                     plot ([hLocsRight(1), hLocsRight(1)], [lHeight - lWidth, 0], '--', 'color', [0.4, 0.4, 0.4]);
%                                     plot ([hLocsLeft(end), hLocsLeft(end)], [lHeight - lWidth, 0], '--', 'color', [0.4, 0.4, 0.4]);

                                end
                            end
                        end
                    else
                        scatter(hLocsLeft(end) + 1, 0, 40, 'b', 'Marker', 'o', 'Linewidth', 2)

                        if numFromBound == 1
                            scatter(hLocsRight(1) - 1, 0, 400, 'r', 'Marker', '.', 'Linewidth', 2)
                        end
                        plot([hLocsLeft(end), hLocsLeft(end) + h * N], [0, 0], 'b--', 'Linewidth', 2)
    %                     xtickLocs = [hLocsLeft(end), hLocsLeft(end) + h * N, hLocsRight(1)];
    %                     xlabelSave = ["$u_M$", "$I_3'\mathbf{v}^n$", "$w_0$"];        
                        if hLocsLeft(end) == hLocsRight(1)
                            xtickLocs = [hLocsRight(1) - 1, hLocsLeft(end), hLocsLeft(end) + 1, hLocsRight(1) + 2];
                            xlabelSave = ["$w_{-1}, u_{M-1}$", "$u_M, w_0$", "$\ u_{M+1}, w_1$","$w_2$"]; 
                        else 
                            xtickLocs = [hLocsLeft(end-2), hLocsLeft(end-1), hLocsLeft(end), hLocsRight(1), hLocsRight(2), hLocsRight(3)];
                            xlabelSave = ["$u_{M-2}$", "$u_{M-1}$", "$u_M$", "$w_0$", "$w_1$", "$w_2$"]; 
                         end
                         text(hLocsLeft(end) + h * N, 0.125, "$\ \ I_3\mathbf{v}^n$", 'horizontalalignment', 'center', 'interpreter', 'latex', 'Fontsize', 18);
                        grid on;
                        yticks([0])
    %                     text((hLocsRight(1) + hLocsLeft(end)) * 0.5, 0.1, "$\alpha = "+ num2str(round(alf * 100) / 100) + "$", ...
    %                         'interpreter', 'latex', 'Fontsize', 16, ...
    %                         'horizontalAlignment', 'center', ...
    %                         'color', [0.5, 0.5, 0.5]);
                    end 
                end
            else
                xtickLocs = [0, floor(length(u) / 2), length(u), length(u) + floor(length(w) / 2), N];
                xlabelSave = ["$u_0$", "$u_l$", "$u_M, w_0$", "$w_l$", "$w_{M_w}$"];        
            end
            ylim([-0.6, 0.6])
    %         ylim([-5,5])
            yticks([0])

            xlabel("$l$", 'interpreter', 'latex')
            grid on;
            if zoomed
    %             xlim([hLocsLeft(end-3), hLocsRight(4)])
                if numFromBound == -1
                    xlim([hLocsLeft(13), hLocsRight(4)])
%                     xlim([0.3 * N, 0.7 * N])
%                     xlim([hLocsLeft(end-2) - 0.5 * h * N, hLocsRight(3) + 0.5 * h * N])
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
%             set(gca, 'Fontsize', 16, 'Linewidth', 2,...
%                 'TickLabelInterpreter', 'latex', ...
%                 'xticklabel', xlabelSave, 'XTick', xtickLocs, ...
%                 'Position', [0.02 0.184615384615385 0.950000000000001 0.801030548398969]);
% %                 'Position', [0.02 0.129186602870813 0.950000000000001 0.85645933014354]);
            set(gcf, 'Color', 'w');

            if n == 16
               disp("wait")
            end
        else
            
            if plotROCEnergy
                subplot(311)

                plot([1:length(u), (1:length(w))+length(u)-1], [u; w])
                subplot(312)
                plot(totROCenergy(1:n));

                subplot(313)
                hold off;
                plot (rOCconnectionU(1:n));
    %             plot(rOCKinEnergyU(1:n));
                hold on;
                plot (rOCconnectionW(1:n));

    %             plot(rOCPotEnergyU(1:n));
    %             plot(rOCKinEnergyW(1:n));
    %             plot(rOCPotEnergyW(1:n));
                pause(0.1)
            else
                subplot(311)
                plot([1:length(u), (1:length(w))+length(u)-1], [u; w])

                subplot(312)
                plot(totTotEnergy(1:n) / totTotEnergy(1) - 1);
%                 plot(totTotEnergy(1:n))
                subplot(313)
                hold off;
                plot((totEnergy(2:n) - totEnergy(1)));

                hold on;
                plot(-boundaryEnergyUH(1:n-1)-boundaryEnergyWH(1:n-1));

            end
        drawnow;
%         if frame <= length(Mvid) && filmFlag == true
%             Mvid(frame) = getframe(gcf);
%             frame = frame + 1;
%         else
%             v = VideoWriter('dynamicGridCenterZoomedOut.mp4', 'MPEG-4');
%             v.FrameRate = 15;
%             open(v)
%             writeVideo(v, Mvid);
%             close(v)
%             filmFlag = false;
%             break;
%         end
        end
        
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
end
if ~setting
    figure
    spectrogram(outputD,512,64,512, curFs, 'yaxis');
%     set(gcf, 'color', 'w')
end
% subplot(2,1,1)
% hold on
% plot((1:lengthSound) / fs, outFree)

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