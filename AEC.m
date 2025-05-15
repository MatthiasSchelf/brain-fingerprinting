function AEC = AEC(a)
% computes Amplitude Envelope Correlation with pairwise orthigonalization
% a is a filtered multichannel signal (time x channels)
% Optimized implementation with vectorization and warning suppression

% Get dimensions
[T, N] = size(a);
AEC = zeros(N, N);

% Precompute envelope for all channels at once
amp = zeros(size(a));
for ch = 1:N
    amp(:, ch) = envelope(a(:, ch), 1, 'analytic');
end

% Vectorized computation for channel pairs
% Preallocate arrays for vectorized operations
ort_signals = zeros(T, N, N);  % Store orthogonalized signals
amp_ort = zeros(T, N, N);      % Store their envelopes

% Precompute all orthogonalizations (can be done outside the loops)
for i = 1:N
    for j = 1:N
        if i ~= j
            ort_signals(:, i, j) = orthog_timedomain(a(:, i), a(:, j));
            amp_ort(:, i, j) = envelope(ort_signals(:, i, j), 1, 'analytic');
        end
    end
end

% Calculate correlations and fill the AEC matrix
for i = 1:N
    for j = i+1:N
        % Get correlation in both directions
        corr1 = corrcoef(amp_ort(:, i, j), amp(:, i));
        corr2 = corrcoef(amp_ort(:, j, i), amp(:, j));

        % Take absolute values and average
        AEC_mean = (abs(corr1(1,2)) + abs(corr2(1,2))) / 2;

        % Fill symmetric matrix
        AEC(i, j) = AEC_mean;
        AEC(j, i) = AEC_mean;
    end
end

end

function [Yorth]=orthog_timedomain(X,Y)
R=[ones(length(X),1) X] \ Y;
Ypred=X*R(2);
Yorth=Y-Ypred;

end