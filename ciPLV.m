function ciPLV = ciPLV(filtered)
% Computes the corrected imaginary Phase Locking Value (ciPLV)
% filtered: Input signals (time x channels)

% Apply Hilbert transform to get analytic signals
% Note: Hilbert operates on columns, so we transpose after to get channels x time
analytic_signals = hilbert(filtered)';

% Get number of time points for normalization
[~, num_timepoints] = size(analytic_signals);

% Normalize to get phase (unit complex numbers on the unit circle)
phase = analytic_signals ./ abs(analytic_signals);

% Compute cross-spectral density matrix
csd = phase * phase';

% Pre-compute the normalization factor (1/n)
inv_n = 1 / num_timepoints;

% Calculate the normalized real and imaginary parts
real_csd_norm = real(csd) * inv_n;
imag_csd_norm = imag(csd) * inv_n;

% Calculate ciPLV using the formula: |Im(CSD)/n| / sqrt(1-(Re(CSD)/n)^2)
ciPLV = abs(imag_csd_norm ./ sqrt(1 - real_csd_norm.^2));
end