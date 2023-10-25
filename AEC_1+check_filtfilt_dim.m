%%amplitude_envelope_correlation
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

%% Make a plot to check if the filtering is done correctly.

sensor = 1; % Can be any sensor really. 

% Create a figure
figure;

% Plot the original data in blue
subplot(2, 1, 1);
plot(1:6000, sensors_sub01_ses01(sensor, :), 'b');
title(['Original Data - Sensor ', num2str(sensor)]);
xlabel('Time');
ylabel('Activation');

% Plot the filtered data in red
subplot(2, 1, 2);
plot(1:6000, beta_sub01ses01(:, :), 'r');
title(['Filtered Data - Sensor ', num2str(sensor)]);
xlabel('Time');
ylabel('Activation');

% This looks good. This is, what I would assume the right way to go. 

%%

%Let's now transpose the matrices so that time is on the first dimension.

sensors_sub01_ses01 = sensors_sub01_ses01.';
sensors_sub01_ses02 = sensors_sub01_ses02.';
sensors_sub02_ses01 = sensors_sub02_ses01.';
sensors_sub02_ses02 = sensors_sub02_ses02.';
sensors_sub03_ses01 = sensors_sub03_ses01.';
sensors_sub03_ses02 = sensors_sub03_ses02.';

%%

% Let's filter the data again. 

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

%% Make a plot to check if the filtering is done correctly. 

%Comparing the raw and filtered data. 

sensor = 1; %Can be any sensor really. 

% Create a figure
figure;

% Plot the original data in blu
subplot(2, 1, 1);
plot(1:273, sensors_sub01_ses01(sensor, :), 'b'); % Here we put in 273 instead of
% 6000 because otherwise the code won't work.
title(['Original Data - Sensor ', num2str(sensor)]);
xlabel('Time');
ylabel('Activation');

% Plot the filtered data in red
subplot(2, 1, 2);
plot(1:273, beta_sub01ses01(:, :), 'r');% Here we put in 273 instead of
% 6000 because otherwise the code won't work.
title(['Filtered Data - Sensor ', num2str(sensor)]);
xlabel('Time');
ylabel('Activation');

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

