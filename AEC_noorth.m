function AEC = AEC_noorth(a)
% a is a filtered multichannel signal (time x channels)
% Optimized implementation with vectorization

% Get dimensions
[~, N] = size(a);

% Vectorized computation of amplitude envelopes
% (We'll keep the loop for envelope computation as it may not be vectorizable
% depending on how envelope() is implemented)
amp = zeros(size(a));
for ch = 1:N
    amp(:, ch) = envelope(a(:, ch), 1, 'analytic');
end

% Suppress rank deficiency warnings
warning('off', 'MATLAB:rankDeficientMatrix');

% Compute all correlations at once
% corrcoef(X) returns correlation matrix between all columns of X
corr_matrix = corrcoef(amp);

% Take absolute values of all correlations
AEC = abs(corr_matrix);

% Restore warning state
warning('on', 'MATLAB:rankDeficientMatrix');
end