close all;
firstIteration = true;

if firstIteration
    clear all;
    firstIteration = true;
else
    figure('Position', [200, 200, 1200, 450])
    set(gcf, 'color', 'w');
end
% figure;
loopingN = false;
loopNStart = 50.0; % also use for Ninit an Nend
loopNend = 60.0;
plotModeShapesBool = ~firstIteration;

lowPassConnection = false;
lpExponent = 30;

modeToPlot = 8; % if -1 plot all modes
limSubplots = 3;
if loopingN
    range = loopNStart:0.01:loopNend;
else
    range = 1;
end
if plotModeShapesBool
    loopAmount = 100;
elseif loopingN
    loopAmount = 1000;
else
    loopAmount = 10000;
end
if firstIteration
    loopAmountRange = 1:loopAmount;
else
    loopAmountRange = 1:(loopAmount+1);
end
interpolation = "sinc";

%{
    0: fullSinc is false 
    1: include all moving points
    2: include boundaries as well
    3: include virtual grid points
%}
fullSinc = 3; 

%{
    Number from the right boundary (quite important, switches between
    different techniques)
    -1: Adding to the center alternating between left and right string.
    0: Interpolated boundary
    1: Right string has a single moving point. Using simply supported boundary condition
    2: Right string has two moving points. When trying to solve the cubic
    interpolation, w_2 is always 0 (that's why this is a bit different)
    >3: (Expected behaviour) Selects where to add points (to left string).
%}

numFromBound = -1;

for Nloop = range
    fs = 44100;                 % sample rate
    k = 1/fs;                   % time step
    
    if loopingN
        Ninit = Nloop;          % number of intervals (!)
        Nend = Nloop;
    else
        Ninit = loopNStart;
        Nend = loopNend;
    end
    
    h = 1/Ninit;                % grid spacing
    c = h/k;                    % wave speed
    N = floor(1/h);             % recalculate (should be exactly equal to Ninit)
    NinitSave = Ninit;
    
    cInit = h/k;
    cEnd = 1/(Nend * k);
    cVec = linspace (cInit, cEnd, loopAmount + 1);
    
    if interpolation == "none"
        h = 1/N;
    end
    
    lambdaSq = (cInit * k / h)^2;
    
    
    % Create B-matrix
    BFull = zeros(N, N); % matrix for scheme without boundaries, plus one overlap
    
    
    if numFromBound == -1 % add alternatively to both in the center
        M = ceil(N/2);
        Mw = floor(N/2);
    else
        M = N-numFromBound;
        Mw = numFromBound;
    end
    
    if interpolation == "none"
        Dxxu = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        BFull(1:N, 1:N) = 2 * eye(N) + lambdaSq * Dxxu;
    else
        Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
            sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
            sparse(1:M-1, 2:M, ones(1, M-1), M, M));
        BFull(1:M, 1:M) = 2 * eye(M) + lambdaSq * Dxxu;

        Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
            sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
            sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));

        BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + lambdaSq * Dxxw;

        BFullInit = BFull;
    end
    maxNumberOfPoints = max(Ninit, Nend);
    if firstIteration
        modesSave = zeros(loopAmount, floor(maxNumberOfPoints));
    end
    NPrev = N;
    j = 1;
    if firstIteration
        loopStart = zeros(ceil(abs(Nend-Ninit)), 1);
        loopStart(j) = 1;
    end
    addLastPoint = true;
    
    for i = loopAmountRange
        h = cVec(i) * k;
        Ninit = 1/h;
        N = floor(Ninit);
        alf  = Ninit - N;
        alfSave(i) = alf;
        if interpolation == "none"
            h = 1/N;
        end
        
        lambdaSq = (cVec(i) * k / h)^2;
        
        if N ~= NPrev
%             modeToPlot = modeToPlot + 1;
            BFull = zeros(N, N);
            if i~= 1 && firstIteration
                j = j + 1;
                loopStart(j) = i;
            end
            if i == loopAmount
                addLastPoint = false;
            end
            
            if numFromBound == -1 % add alternatively to both in the center
                M = ceil(N/2);
                Mw = floor(N/2);
            else
                M = N-numFromBound;
                Mw = numFromBound;
            end
            
            Dxxu = (sparse(2:M, 1:M-1, ones(1, M-1), M, M) + ...
                sparse(1:M, 1:M, -2 * ones(1, M), M, M) + ...
                sparse(1:M-1, 2:M, ones(1, M-1), M, M));
            BFull(1:M, 1:M) = 2 * eye(M) + lambdaSq * Dxxu;
            
            Dxxw = (sparse(2:Mw, 1:Mw-1, ones(1, Mw-1), Mw, Mw) + ...
                sparse(1:Mw, 1:Mw, -2 * ones(1, Mw), Mw, Mw) + ...
                sparse(1:Mw-1, 2:Mw, ones(1, Mw-1), Mw, Mw));
            
            BFull((M+1):end, (M+1):end) = 2 * eye(Mw) + lambdaSq * Dxxw;
            
            BFullInit = BFull;
        end
        NPrev = N;
        
        if interpolation == "none"
            
            Dxxu = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
            BFull(1:N, 1:N) = 2 * eye(N) + lambdaSq * Dxxu;
            
        elseif interpolation == "linear"
            BFull(M, (M+1):(length(u) + 2)) = [alf, (1-alf)];
            BFull(M + 1, (M-1) : M) = [(1-alf), alf];
        elseif interpolation == "cubic"
            ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
            Ainv = inv([1, -ip(4); -ip(4), 1]);
            if numFromBound == 2
                BFull(M, (M + 1):(M + 2)) = ip(3:-1:2) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 2)) = BFullInit(M+1, (M + 1):(M + 2)) + ip(3:-1:2) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            elseif numFromBound == 1
                BFull(M, M+1) = (ip(3) - ip(1)) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, M+1) = BFullInit(M+1, M+1) + (ip(3) - ip(1)) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            elseif numFromBound == 0
                hLocs = 1:-h:0;
                if hLocs(end) >= h/2
                    alf = (h-hLocs(end))/ hLocs(end);
                    %                     uVirtual = -alf * u(1);
                    BFull(N, N) = BFullInit(N, N) - alf;
                    
                else
                    alf = (2*hLocs(end)) / h;
                    %                     uVirtual = -(alf * u(1) + (1-alf) * u(2));
                    BFull(N, N-1:N) = BFullInit(N, N-1:N) - [(1-alf), alf];
                end
                %                 BFull(N, N) = -1;
            else
                BFull(M, (M + 1):(M + 3)) = ip(3:-1:1) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + ip(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 3)) = BFullInit(M+1, (M + 1):(M + 3)) + ip(3:-1:1) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            end
            
        elseif interpolation == "sinc"
            includeUMp1AndWm1 = true;
            if alf < 1e-6
                alf = alf + 1e-6;
            end
            
            alphaBand = 1; % relative bandwidth range
            bmax = alphaBand*pi;
            
%             sincWidth = 2;
            sincWidth = floor(N / 2) - 1;
%             sincWidth = numFromBound+1;
            if fullSinc == 0
                if includeUMp1AndWm1
                    xUMp1 = [-sincWidth:-1, -1:sincWidth-1]';
                    xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth+1, 1)];
                else
                    xUMp1 = (-sincWidth:sincWidth-1)';
                    xUMp1 = xUMp1 + [zeros(sincWidth, 1); alf * ones(sincWidth, 1)];
                end
            elseif fullSinc == 2
                xUMp1 = (0:N+1)' - M - 1;
                xUMp1 = xUMp1 + [zeros(M+1, 1); alf * ones(N-M+1, 1) - 1];
            elseif fullSinc == 3
                xUMp1 = (-1:N+2)' - M - 1;
                xUMp1 = xUMp1 + [zeros(M+2, 1); alf * ones(N-M+2, 1) - 1];
%                 xUMp1 = (-1:N+2)' - M - 2;
%                 xUMp1 = xUMp1 + [zeros(M+3, 1); alf * ones(N-M+1, 1) - 1];
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
%             if alf == 0
%                 idxAu = find(round(aU) == 1);
%                 aU(idxAu-2) = -1;
%                 aU(idxAu-1) = 1;
%                 aU(idxAu) = 1;
%             end
            if fullSinc == 0
                if includeUMp1AndWm1
                    xWm1 = [-sincWidth+1:1, 1:sincWidth]';
                    xWm1 = xWm1 - [alf * ones(sincWidth+1, 1); zeros(sincWidth, 1)];
                else
                    xWm1 = (-sincWidth+1:sincWidth)';
                    xWm1 = xWm1 - [alf * ones(sincWidth, 1); zeros(sincWidth, 1)];
                end
            elseif fullSinc == 3
                xWm1 = (-1:N+2)' - M;
                xWm1 = xWm1 - [alf * ones(M+2, 1) - 1; zeros(N-M+2, 1)];
%                 xWm1 = (-1:N+2)' - M;
%                 xWm1 = xWm1 - [alf * ones(M+3, 1); ones(N-M+1, 1)];
            end
            bW = (sin(bmax*xWm1)./xWm1);
            if sum(isnan(bW))
                bW(isnan(bW)) = bmax;
            end 
            distW = xWm1*ones(1,iLen)-ones(iLen,1)*xWm1';    % distance matrix between points
            AW = sin(bmax*distW)./distW;
            AW(1+(iLen+1)*[0:iLen-1]') = bmax;         % collection of sinc functions with centers at grid point locations
            AW(isnan(AW)) = bmax;
            aW = AW\bW; %optimal coefficients

            if fullSinc == 3
                BFull(M, :) = BFullInit(M, :) + aU(3:end-2)' - [aU(1), zeros(1, length(aU)-6), aU(end)];
                BFull(M+1, :) = BFullInit(M+1, :) + aW(3:end-2)' - [aW(1), zeros(1, length(aW)-6), aW(end)];
            elseif M + sincWidth == N+1 % use boundary condition
                if includeUMp1AndWm1
                    BFull(M, M-sincWidth+1:end) = BFullInit(M, M-sincWidth+1:end) + aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)];
                    BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW(1:end-1)';
                else
                    BFull(M, M-sincWidth+2:end) = BFullInit(M, M-sincWidth+2:end) + aU(1:end-2)' - [zeros(1, length(aU)-3), aU(end)];
                    BFull(M+1, (M-sincWidth: end-1)) = BFullInit(M+1, (M-sincWidth: end-1)) + aW(1:end-1)';
                end
            elseif M + sincWidth == N % use boundary condition
                if includeUMp1AndWm1
                    BFull(M, M-sincWidth+1:end) = BFullInit(M, M-sincWidth+1:end) + aU(1:end-1)';
                    BFull(M+1, (M-sincWidth : end)) = BFullInit(M+1, (M-sincWidth : end)) + aW';
                else
                    BFull(M, M-sincWidth+2:end) = BFullInit(M, M-sincWidth+2:end) + aU(1:end-1)';
                    BFull(M+1, (M-sincWidth: end-1)) = BFullInit(M+1, (M-sincWidth: end-1)) + aW';
                end
            else
%                 BFull(M, (M + 1):(M + 3)) = aU(1:3);
                if includeUMp1AndWm1
                    sincRange = (M-sincWidth : M+sincWidth);
                    BFull(M, sincRange+1) =  BFullInit(M, sincRange+1) + aU';
                    BFull(M+1, sincRange) = BFullInit(M+1, sincRange) + aW';
                else
                    sincRange =[M-sincWidth:M-1, M+1:M+sincWidth];
                    BFull(M, sincRange+1) =  BFullInit(M, sincRange+1) + aU';
                    BFull(M+1, sincRange) = BFullInit(M+1, sincRange) + aW';
                end
%                 if  mod(i, 100) == 0 && i>0 
%                     hold off   
%                     plot(aU)
%                     hold on;
%                     plot(aW)
% %                     pause(0.5)
%                     title(alf)
%                     ylim([-1,1])
%                     drawnow;
%                 end
%                 BFull(M+1, (M-2 : M)) =  ip(1:3) * Ainv(2, 2);
            end
            if lowPassConnection
                lpVec = 0.5 * [-(1-alf)^(lpExponent), (1-alf)^(lpExponent)];
                BFull(M, M:M+1) = BFull(M, M:M+1) + lpVec;
                BFull(M+1, M:M+1) = BFull(M+1, M:M+1) - lpVec;
            end
%             if mod(i, 1) == 0
%                 subplot(2,1,1)
%                 plot (aU)
%                 subplot(2,1,2)
%                 plot (aW)
%                 drawnow;
%             end
        elseif interpolation == "altSinc"
            alphaBand = 1; % relative bandwidth range
            bmax = alphaBand*pi / h;
            
            xIpLocs = (-alf + (-1:2)) * h;
            sincIp = sin(bmax * xIpLocs) ./ (bmax * xIpLocs);
            sincIp(isnan(sincIp)) = 1;
%             sincIp = fliplr(sincIp);
            Ainv = inv([1, -sincIp(4); -sincIp(4), 1]);
            
            if numFromBound == 1
                BFull(M, M+1) = (sincIp(3) - sincIp(1)) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + sincIp(1:3) * Ainv(1, 2);
                BFull(M+1, M+1) = BFullInit(M+1, M+1) + (sincIp(3) - sincIp(1)) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  sincIp(1:3) * Ainv(2, 2);
                
            else
                
                BFull(M, (M + 1):(M + 3)) = sincIp(3:-1:1) * Ainv(1, 1);
                BFull(M, (M-2 : M)) =  BFullInit(M, (M-2) : M) + sincIp(1:3) * Ainv(1, 2);
                BFull(M+1, (M + 1):(M + 3)) = BFullInit(M+1, (M + 1):(M + 3)) + sincIp(3:-1:1) * Ainv(2, 1);
                BFull(M+1, (M-2 : M)) =  sincIp(1:3) * Ainv(2, 2);
            end            
            if mod(i, 100) == 0
                hold off;
                plot (xIpLocs, sincIp)
                hold on;
                plotRange = (min(xIpLocs) : 0.001 : max(xIpLocs));
                plot (plotRange, sin(bmax * plotRange) ./ (bmax * plotRange))
                title(alf)
                drawnow;
            end
        end
        % imagesc(BFull)
        % drawnow;          
        [~, D, W] = eig(BFull, 'vector');
        if plotModeShapesBool
    %         order = Ninit:-1:1;
            plotModeShapes;
        end
        if firstIteration
            modesSave(i, 1:N) = sort(1/(2 * pi * k) * acos (1/2 * D));
        end
    end
    if addLastPoint && firstIteration
        loopStart(j+1) = loopAmount;
    end
    
    if firstIteration
        hold off;
        plotModesSave;
    end

    drawnow;
end
xData = 1:loopAmount;
goodness = zeros(floor(loopNStart), 1);
rms = zeros(floor(loopNStart), 1);
sse = zeros(floor(loopNStart), 1);

% figure;
% for i = 1:loopNStart
%     [~, got] = fit (xData(~isnan(modesSave(:,i)))', modesSave(~isnan(modesSave(:,i)),i), 'poly1');
%     goodness(i) = got.adjrsquare;
%     rms(i) = got.rmse;
%     sse(i) = got.sse;
%     
%     %     got
%     %     eval(['mode', num2str(i), ' = modesSave(:,', num2str(i), ');']);
% end
% subplot(211)
% plot(goodness)
% subplot(212)
% plot(sse)


