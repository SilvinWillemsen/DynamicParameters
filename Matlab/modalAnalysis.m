% close all;
clear all;
figure;
loopingN = false;
loopNStart = 15.0; % also use for Ninit an Nend
loopNend = 20;

if loopingN
    range = loopNStart:0.01:loopNend;
else
    range = 1;
end
loopAmount = 10000;
interpolation = "none"; %!!!!

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
    
    lambdaSq = (cInit * k / h)^2
    
    
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
    modesSave = zeros(loopAmount, floor(maxNumberOfPoints));
    NPrev = N;
    j = 1;
    loopStart = zeros(ceil(abs(Nend-Ninit)), 1);
    loopStart(j) = 1;
    addLastPoint = true;
    
    for i = 1:loopAmount
        h = cVec(i) * k;
        Ninit = 1/h;
        N = floor(Ninit);
        alf  = Ninit - N;
        
        if interpolation == "none"
            h = 1/N;
        end
        
        lambdaSq = (cVec(i) * k / h)^2;
        
        if N ~= NPrev
            BFull = zeros(N, N);
            if i~= 1
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
        else
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
            
        end
        % imagesc(BFull)
        % drawnow;
        modesSave(i, 1:N) = sort(1/(2 * pi * k) * acos (1/2 * eig(BFull)));
    end
    if addLastPoint
        loopStart(j+1) = loopAmount;
    end
    %     figure
    hold off;
    modesSave = modesSave;
    modesSave(modesSave==0) = nan;
    h = plot(real(modesSave));
    
    colours = [];
    for j = 1:floor(maxNumberOfPoints)
        if mod(j,2) == 0
            colours = [colours; 0,0,1];
        else
            colours = [colours; 1,0,0];
        end
    end
    %     figure
    set(h, {'color', 'Linewidth'}, [num2cell(colours, 2), num2cell(2 * ones(floor(maxNumberOfPoints), 1))])
    title ("Modal Analysis $N = " + loopNStart + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
    xlabelsave = num2cell(NinitSave:Nend);
    set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStart, 'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex')
    xlabel("$N$", 'interpreter', 'latex')
    ylabel("Frequency (Hz)", 'interpreter', 'latex')
    ylim([0, fs / 2])
    grid on
    drawnow;
end
xData = 1:loopAmount;
goodness = zeros(loopNStart, 1);
rms = zeros(loopNStart, 1);
sse = zeros(loopNStart, 1);

figure;
for i = 1:loopNStart
    [~, got] = fit (xData(~isnan(modesSave(:,i)))', modesSave(~isnan(modesSave(:,i)),i), 'poly1');
    goodness(i) = got.adjrsquare;
    rms(i) = got.rmse;
    sse(i) = got.sse;
    
    %     got
    %     eval(['mode', num2str(i), ' = modesSave(:,', num2str(i), ');']);
end
subplot(211)
plot(goodness)
subplot(212)
plot(sse)


