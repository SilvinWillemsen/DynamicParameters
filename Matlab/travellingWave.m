clear all;

fs = 44100;
k = 1/fs;
lengthSound = fs;
N = 30;           % edit how many points you want
h = 1/N;

cInit = h/k;            % calculate wave speed
c = cInit;

% uLeft = zeros(N+1, 1);
% uRight = zeros(N+1, 1);
% uLeft(floor(N/5)-4:floor(N/5)+4) = 0.5 * hann(9);
% uRight(floor(N/5)-4:floor(N/5)+4) = 0.5 * hann(9);

width = 10;
excitation = (1 - cos(2 * pi * [0:width] / width)) * 0.5;

for n = 1:lengthSound
    hold off;
    
%     plot (uLeft, 'Linewidth', 2)
%     
%     hold on;
%     plot (uRight, 'Linewidth', 2)
    plot(uLeft + uRight, 'Linewidth', 2)
    
    uLeftBound = uLeft(1);
    uRightBound = uRight(end);
    uLeft = [uLeft(2:end); -uRightBound];
    uRight = [-uLeftBound; uRight(1:end-1)];

    pause(0.1);
    drawnow;
end
