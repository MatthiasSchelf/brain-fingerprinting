%% amplitude_envelope_correlation
% of given sensordata, head models and atlases


%% Filter the data so only beta frequencies remain (like the paper of Da Silva et al. (2021))
% (da Silva Castanheira, J., Orozco Perez, H.D., Misic, B. et al. Brief
% segments of neurophysiological activity enable individual differentiation. Nat Commun 12,
% 5713 (2021). https://doi.org/10.1038/s41467-021-25895-8)

%%

%Start by transposing the matrices so that time is on the first dimension.

sensors_sub01_ses01 = sensors_sub01_ses01.';
sensors_sub01_ses02 = sensors_sub01_ses02.';
sensors_sub02_ses01 = sensors_sub02_ses01.';
sensors_sub02_ses02 = sensors_sub02_ses02.';
sensors_sub03_ses01 = sensors_sub03_ses01.';
sensors_sub03_ses02 = sensors_sub03_ses02.';

%% Filter the data so only beta frequencies remain (like the paper of Da Silva et al. (2021))
% (da Silva Castanheira, J., Orozco Perez, H.D., Misic, B. et al. Brief
% segments of neurophysiological activity enable individual differentiation. Nat Commun 12,
% 5713 (2021). https://doi.org/10.1038/s41467-021-25895-8)

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

%Define a one second segment
start_time = 3000;
end_time = 3600;

% Define the same axis limits for all subplots
axis_limits = [1, 6000, min(beta_sub01ses01(:, 1)), max(beta_sub01ses01(:, 1))];

% Create a figure
figure;

% Plot the whole channel
subplot(1,3,1);
plot(1:6000, beta_sub01ses01(1:6000, 1), 'b'); % Use ':' to include all time points
title(['Time Course of Channel ', num2str(1)]);
xlabel('Time');
ylabel('Activation');
axis(axis_limits);

% Plot the whole channel excluding the first and last 100 time points. 
subplot(1, 3, 2);
plot(101:5900, beta_sub01ses01(101:5900, 1), 'b');
title('Time Course - end and beginning');
xlabel('Time');
ylabel('Activation');
axis(axis_limits);

%Plot only a one second segment
subplot(1, 3, 3);
plot(start_time:end_time, beta_sub01ses01(start_time:end_time, 1), 'b');
title(['1-Second Segment of Channel ', num2str(1)]);
xlabel('Time');
ylabel('Activation');

%% Let's check by putting all sensors to zero except one.

% Create a new matrix where all channels except channel 1 are set to 0
check_matrix = sensors_sub01_ses01;
check_matrix(:, 2:end) = 0; % Set all columns except the first one to 0

%% Filter

beta_check = filtfilt(b, a, check_matrix);

%% Now plot the matrix in which only one channel is present
%Define a one second segment
start_time = 3000;
end_time = 3600;

% Define the same axis limits for all subplots
axis_limits = [1, 6000, min(beta_check(:, 1)), max(beta_check(:, 1))];

% Create a figure
figure;

% Plot the whole channel
subplot(1,3,1);
plot(1:6000, beta_check(1:6000), 'b'); % Use ':' to include all time points
title(['Time Course of Channel - Check ', num2str(1)]);
xlabel('Time');
ylabel('Activation');
axis(axis_limits);

% Plot the whole channel excluding the first and last 100 time points. 
subplot(1, 3, 2);
plot(101:5900, beta_check(101:5900), 'b');
title('Time Course - end and beginning - check');
xlabel('Time');
ylabel('Activation');
axis(axis_limits);

%Plot only a one second segment
subplot(1, 3, 3);
plot(start_time:end_time, beta_check(start_time:end_time), 'b');
title(['1-Second Segment of Channel - check', num2str(1)]);
xlabel('Time');
ylabel('Activation');

%% Combining sensordata and head model + transposing the matrix again. 
 
sub01ses01head01 = kernel_sub01 * beta_sub01ses01.';

sub01ses02head01 = kernel_sub01 * beta_sub01ses02.';

sub01ses01head02 = kernel_sub02 * beta_sub01ses01.';

sub01ses02head02 = kernel_sub02 * beta_sub01ses02.';

sub02ses01head02 = kernel_sub02 * beta_sub02ses01.';

sub02ses01head01 = kernel_sub01 * beta_sub02ses01.';

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


