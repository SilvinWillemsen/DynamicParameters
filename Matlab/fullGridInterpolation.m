%{
    Full Grid Interpolation
    
    Function output is a waveform
    Function arguments
    - N start for fs = 44100;
    - N end for fs = 44100;
    - sample rate
    - output length (in s)
    - excitation location center as ratio of length 
    - excitation width as ratio of length
    - location of the output as a ratio of length
%}

% function [out] = fullGridInterpolation (startN, endN, fs, outLength, excitationLoc, excitationWidth, outputLocStart)
    

clear all;
close all;
clc; 

load ("/Users/SilvinW/Desktop/output/FullGridInterpolation-ejjbtghhbugjpfgqiqqcdpfkqwvo/Build/Products/Debug/stateAt.csv")
hold off; plot((1:29) / 30, stateAt(1:29));
hold on;
plot((1:30) / 31, stateAt(30:end));


startN = 30;
endN = 30;
fs = 44100;
outLength = 0.5;
excitationLoc = 0.25;
excitationWidth = 0.2;
outputLocStart = 0.9;

lengthSound = outLength * fs;

k = 1/fs;

startC = 44100 / startN;
c = startC;
endC = 44100 / endN;

cDiff = startC - endC;
h = c*k;
N = floor(1/h);
h = 1/N;
NPrev = N;
lambdaSq = c^2 * k^2 / h^2;

    
uNext = zeros(N-1, 1);
u = zeros(N-1, 1);

loc = excitationLoc * N;
width = excitationWidth * N;
raisedCosStart = floor(loc - width  / 2);
raisedCosEnd = floor(loc + width / 2);
u(raisedCosStart:raisedCosEnd) = 0.5 * (1 - cos(2 * pi * [0:width] / width));
uPrev = u;

out = zeros(lengthSound, 1);


outLoc = outputLocStart * N;
NChangeSave = zeros(lengthSound, 1);
idx = 1;

N = N - 1;
Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
    sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
    sparse(1:N-1, 2:N, ones(1, N-1), N, N));
B = 2 * eye(N) + lambdaSq * Dxx;
C = -eye(N);
N = N + 1;
        
for n = 1:lengthSound

    c = startC - n/(lengthSound-1) * cDiff;
    h = c*k;
    N = floor(1/h);
    h = 1/N;
    lambdaSq = c^2 * k^2 / h^2;

    if N ~= NPrev
        if N > NPrev
            u = interp1((1:NPrev-1) / NPrev, u, (1:N-1) / N, 'pchip', 'extrap')';
            uPrev = interp1((1:NPrev-1) / NPrev, uPrev, (1:N-1) / N, 'pchip', 'extrap')';
            uNext = zeros(length(u), 1);

        else
            u = interp1((1:NPrev-1) / NPrev, u, (1:N-1) / N, 'pchip')';
            uPrev = interp1((1:NPrev-1) / NPrev, uPrev, (1:N-1) / N, 'pchip')';
            uNext = zeros(length(u), 1);
        end
        NChangeSave(idx) = n;
        idx = idx + 1;
        
        N = N - 1;
        Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
        B = 2 * eye(N) + lambdaSq * Dxx;
        C = -eye(N);
        N = N + 1;

    end

    NPrev = N;

    uNext = B * u + C * uPrev;

    out(n) = uNext(end - outLoc); % don't need to add interpolation here
%     if mod(n,fs / 44100 * 1000) == 0
    if n == 1001
        hold off;
        plot(u);
        ylim([-1, 1])
        title("$N = " + num2str(N) + "\quad" + num2str((n/lengthSound * 100)) + "\%\ \textrm{done}$", 'interpreter', 'latex')
        hold on;
%         plot(stateAt((1:(N-1)) + (n-1) * (N-1)), '--');
        plot(stateAt, '--')
        pause(0.2)
        drawnow;
    end
    uPrev = u;
    u = uNext;
end
% plot(out)

% end