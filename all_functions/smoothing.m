function Y = smoothing(X, win)

if win == 1
    Y = X;
    return;
end

W = round((win - 1) / 2);
[N, C] = size(X);

Y = zeros(N, C);

if win >= N
    for i = 1 : N
        Y(i, :) = mean(X(max(1, i - W) : min(N, i + W), :));
    end
    return;
end

ts = sum(X(1 : W + 1, :));
for i = 1 : W
    Y(i, :) = ts / (W + i);
    ts = ts + X(W + i + 1, :);
end

for i = W + 1 : N - W - 1
    Y(i, :) = ts / win;
    ts = ts - X(i - W, :) + X(i + W + 1, :);
end

for i = N - W : N
    Y(i, :) = ts / (N - i + W + 1);
    ts = ts - X(i - W, :);
end
    
return
