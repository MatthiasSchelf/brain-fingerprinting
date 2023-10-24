%% amplitude_envelope_correlation
% of given sensordata, head models and atlases

%% Filter the data so only beta frequencies remain (like the paper of Da Silva et al. (2021))
% (da Silva Castanheira, J., Orozco Perez, H.D., Misic, B. et al. Brief
% segments of neurophysiological activity enable individual differentiation. Nat Commun 12,
% 5713 (2021). https://doi.org/10.1038/s41467-021-25895-8)

% First parameters
fs = 600;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

% sub01 ses01

beta_sub01ses01 = filtfilt(b, a, sensors_sub01_ses01);

%sub01 ses02

beta_sub01ses02 = filtfilt(b, a, sensors_sub01_ses02);

% sub02 ses01

beta_sub02ses01 = filtfilt(b, a, sensors_sub02_ses01);

% sub02 ses02 

beta_sub02ses02 = filtfilt(b, a, sensors_sub02_ses02);

% sub03 ses01

beta_sub03ses01 = filtfilt(b, a, sensors_sub03_ses01);

% sub03 ses02 

beta_sub03ses02 = filtfilt(b, a, sensors_sub03_ses02);


%% Combining sensordata and head model
 
sub01ses01head01 = kernel_sub01 * beta_sub01ses01;

sub01ses02head01 = kernel_sub01 * beta_sub01ses02;

sub01ses01head02 = kernel_sub02 * beta_sub01ses01;

sub01ses02head02 = kernel_sub02 * beta_sub01ses02;

sub02ses01head02 = kernel_sub02 * beta_sub02ses01;

sub02ses01head01 = kernel_sub01 * beta_sub02ses01;

%% Using the Hilber transform 

envelope_sub01ses01head01 = abs(hilbert(sub01ses01head01));
envelope_sub01ses02head01 = abs(hilbert(sub01ses02head01));
envelope_sub01ses01head02 = abs(hilbert(sub01ses01head02));
envelope_sub01ses02head02 = abs(hilbert(sub01ses02head02));
envelope_sub02ses01head02 = abs(hilbert(sub02ses01head02));
envelope_sub02ses01head01 = abs(hilbert(sub02ses01head01));

%% Compute the AEC. 

% Just a small test, this should be = 1.
aec_value0 = corr2(envelope_sub01ses01head01, envelope_sub01ses01head01);

% Just between sessions. Test-restest bestween two different sessions. 
aec_value1 = corr2(envelope_sub01ses01head01, envelope_sub01ses02head01);

% Only changing the head to see the influence of the head model.
aec_value2 = corr2(envelope_sub01ses01head01, envelope_sub01ses01head02);

% Only changing the sensor data but using the same headmodel. 
aec_value3 = corr2(envelope_sub01ses01head01, envelope_sub02ses01head01);

% Changing head model and session. 
aec_value4 = corr2(envelope_sub01ses01head01, envelope_sub01ses02head02);

