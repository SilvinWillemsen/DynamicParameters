clear all;
close all;
clc;

% drawing variables
drawSpeed = 1000;
drawStart = 1;

fs = 44100;             % Sample rate
k = 1/fs;               % Time step
lengthSound = fs*2;     % Length of the simulation

Ninit = 30.5;           % edit how many points you want
h = 1/Ninit;

cInit = h/k;            % calculate wave speed
c = cInit;

h = c*k;                % calculate h from wavespeed
N = floor(1/h);         % calculate points from h

lambdaSq = (c*k/h)^2    % should always be 1 as h is not recalculated

alf = Ninit - N;        % fractional remainder for the grid point

%% initialise states
uNext = zeros(ceil(N/2), 1);
u = zeros(ceil(N/2), 1);

u(floor(N/5)-4:floor(N/5)+4) = hann (9); % use hann window for excitation
uPrev = u;

wNext = zeros(floor(N/2), 1);
w = zeros(floor(N/2), 1);
wPrev = w;

% initialise laplacians
eu = ones(length(u), 1);
Dxxu = spdiags([eu -2*eu eu], -1:1, length(u),length(u));

ew = ones(length(w), 1);
Dxxw = spdiags([ew -2*ew ew], -1:1, length(w),length(w));

interpol = "cubic";
outFree = zeros(lengthSound, 1);

changeC = true; % set to true for dynamic changes in wavespeed

%% main loop
for n = 1:lengthSound  
 
    % change wave speed...
    if changeC
        c = cInit * (1-0.2 * sin(15 * pi * n/fs));  % very extremely 
    else
        c = c;
    end
    cSave(n) = c;
    
    % save previous state for comparison later
    NPrev = N;

    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambdaSq = c^2 * k^2 / h^2;

    alf = (Ninit - N);
  
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
    interplatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                        ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];
    
    %% left half string
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    uNext(end) = uNext(end) + lambdaSq * interplatedPoints(1);
    
    %% right half string
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + lambdaSq * interplatedPoints(2);

    %% save output
    outFree(n) = w(end - 5);

    %% draw stuff
    if mod(n, drawSpeed) == 0 && n > drawStart
        hLocsLeft = (1:(length(u))) * h;
        hLocsRight = flip(1 - ((1:(length(w))) * h));
   
        hold off;
        plot(hLocsLeft, u, 'LineWidth' ,2, 'Marker', '.', 'MarkerSize', 20, 'Color', 'b') 

        hold on;
        wOffset = 0.00;
        plot(hLocsRight, w + wOffset, 'Linewidth', 2,  'Marker', 'o', 'MarkerSize', 10, 'Color', 'r')
                
        ylim([-1.5, 1.5])
        grid on;
        title("Sample = " + num2str(n) + "   N = " + num2str(floor(Ninit * 10) / 10))
        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        if n == 16
            disp("wait")
        end
        drawnow;
        
    end
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
end

% subplot(2,1,1)
hold off
plot(outFree)