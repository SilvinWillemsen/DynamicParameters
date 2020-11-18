N = 15.0
highQRange = (0:10000) / 10000;
if floor(N) ~= N
    range = ([0:floor(N), N]) / N;
    endLoop = N;

else
    range = (0:N) / N;
    endLoop = N-1;
end

for i = 1:endLoop
    hold off
    plot(range, sin(i * pi * range), '-o');
    hold on;
    plot(highQRange, sin(i * pi * highQRange));

%     drawnow;
%     pause(0.5);
end
N = 10;
Dxx = (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));

[V, D, W] = eig(2 * eye(N) - Dxx);
for i = 1:N
    plot([0; W(:,i); 0])
    drawnow;
    pause(0.5);
end
answer = sort(sqrt(double(solve(test))))

kMat = [1200, -600; -600, 600];
mMat = [2, 0; 0, 2];
detTest = det(kMat - lambdaTest * mMat);
answer = sqrt(double(solve(detTest)))