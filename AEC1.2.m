

%% amplitude_envelope_correlation
% of given sensordata, head models and atlases


%% Filter the data so only beta frequencies remain (like the paper of Da Silva et al. (2021))
% (da Silva Castanheira, J., Orozco Perez, H.D., Misic, B. et al. Brief
% segments of neurophysiological activity enable individual differentiation. Nat Commun 12,
% 5713 (2021). https://doi.org/10.1038/s41467-021-25895-8)

%% Start by transposing the matrices so that time is on the first dimension.

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

%% Here we do the Hilbert transform to get the AEC's from the sensor data.

hil_beta_sub01ses01 = abs(hilbert(beta_sub01ses01));
hil_beta_sub01ses02 = abs(hilbert(beta_sub01ses02));
hil_beta_sub02ses01 = abs(hilbert(beta_sub02ses01));
hil_beta_sub02ses02 = abs(hilbert(beta_sub02ses02));
hil_beta_sub03ses01 = abs(hilbert(beta_sub03ses01));
hil_beta_sub03ses02 = abs(hilbert(beta_sub03ses02));


% Now the data are filtered so that the beta waves are selected and those
% waves are than put through a Hilbert transform to get the AEC

%% Combining sensordata and head model + transposing the matrix again.

% Here we combine the the head model with theAEC data. 
 
sub01ses01head01 = kernel_sub01 * hil_beta_sub01ses01.';

sub01ses02head01 = kernel_sub01 * hil_beta_sub01ses02.';

sub01ses01head02 = kernel_sub02 * hil_beta_sub01ses01.';

sub01ses02head02 = kernel_sub02 * hil_beta_sub01ses02.';

sub02ses01head02 = kernel_sub02 * hil_beta_sub02ses01.';

sub02ses01head01 = kernel_sub01 * hil_beta_sub02ses01.';


%% Combine the sensor within the atlases ROI + plot a connectivity matrix 111

%Create a matrix for the brain activity in the atlas areas
num_areas111 = length(atlas_sub01);
reconstructed111 = zeros(num_areas111, size(sub01ses01head01, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx111 = 1:num_areas111
    area_indices111 = atlas_sub01{area_idx111}; % Get the indices for the current brain area
    area_data111 = mean(sub01ses01head01(area_indices111, :), 1); % Calculate the mean activity within the area
    reconstructed111(area_idx111, :) = area_data111; % Store the area activity in the output matrix
end

% Make connectivity matrices + Hilbert transforms 121
con_matrix111 = corr(reconstructed111');

% Display the correlation matrix
figure;
imagesc(con_matrix111);
colorbar;
title('Connectivity Matrix - Correlation 111');
clim([-1,1]);

%% Combine the sensor within the atlases ROI + plot a connectivity matrix 121

%Create a matrix for the brain activity in the atlas areas
num_areas121 = length(atlas_sub01);
reconstructed121 = zeros(num_areas121, size(sub01ses02head01, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx121 = 1:num_areas121
    area_indices121 = atlas_sub01{area_idx121}; % Get the indices for the current brain area
    area_data121 = mean(sub01ses02head01(area_indices121, :), 1); % Calculate the mean activity within the area
    reconstructed121(area_idx121, :) = area_data121; % Store the area activity in the output matrix
end

% Make connectivity matrix
con_matrix121 = corr(reconstructed121');

% Display the connectivity matrix
figure;
imagesc(con_matrix121);
colorbar;
title('Connectivity Matrix - Correlation 121');
clim([-1,1]);

%% Combine the sensor within the atlases ROI + plot a connectivity matrix 112

%Create a matrix for the brain activity in the atlas areas
num_areas112 = length(atlas_sub02);
reconstructed112 = zeros(num_areas112, size(sub01ses01head02, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx112 = 1:num_areas112
    area_indices112 = atlas_sub02{area_idx112}; % Get the indices for the current brain area
    area_data112 = mean(sub01ses01head02(area_indices112, :), 1); % Calculate the mean activity within the area
    reconstructed112(area_idx112, :) = area_data112; % Store the area activity in the output matrix
end

% Make connectivity matrices + Hilbert transforms 112
con_matrix112 = corr(reconstructed112');

% Display the correlation matrix
figure;
imagesc(con_matrix112);
colorbar;
title('Connectivity Matrix - Correlation 112');
clim([-1,1]);
%% Combine the sensor within the atlases ROI + plot a connectivity matrix 122

%Create a matrix for the brain activity in the atlas areas
num_areas122 = length(atlas_sub02);
reconstructed122 = zeros(num_areas122, size(sub01ses02head02, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx122 = 1:num_areas122
    area_indices122 = atlas_sub02{area_idx122}; % Get the indices for the current brain area
    area_data122 = mean(sub01ses02head02(area_indices122, :), 1); % Calculate the mean activity within the area
    reconstructed122(area_idx122, :) = area_data122; % Store the area activity in the output matrix
end

% Make connectivity matrices + Hilbert transforms 122
con_matrix122 = corr(reconstructed122');

% Display the correlation matrix
figure;
imagesc(con_matrix122);
colorbar;
title('Connectivity Matrix - Correlation 122');
clim([-1,1]);
%% Combine the sensor within the atlases ROI + plot a connectivity matrix 212

%Create a matrix for the brain activity in the atlas areas
num_areas212 = length(atlas_sub02);
reconstructed212 = zeros(num_areas212, size(sub02ses01head02, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx212 = 1:num_areas212
    area_indices212 = atlas_sub02{area_idx212}; % Get the indices for the current brain area
    area_data212 = mean(sub02ses01head02(area_indices212, :), 1); % Calculate the mean activity within the area
    reconstructed212(area_idx212, :) = area_data212; % Store the area activity in the output matrix
end

% Make connectivity matrices + Hilbert transforms 122
con_matrix212 = corr(reconstructed212');

% Display the correlation matrix
figure;
imagesc(con_matrix212);
colorbar;
title('Connectivity Matrix - Correlation 212');
clim([-1,1]);

%% Combine the sensor within the atlases ROI + plot a connectivity matrix 211

%Create a matrix for the brain activity in the atlas areas
num_areas211 = length(atlas_sub01);
reconstructed211 = zeros(num_areas211, size(sub02ses01head01, 2));

%Now we loop through the atlas and the activities to put the right
%activities in the right atlas. 

for area_idx211 = 1:num_areas211
    area_indices211 = atlas_sub01{area_idx211}; % Get the indices for the current brain area
    area_data211 = mean(sub02ses01head01(area_indices211, :), 1); % Calculate the mean activity within the area
    reconstructed211(area_idx211, :) = area_data211; % Store the area activity in the output matrix
end

% Make connectivity matrices + Hilbert transforms 211
con_matrix211 = corr(reconstructed211');

% Display the correlation matrix
figure;
imagesc(con_matrix211);
colorbar;
title('Connectivity Matrix - Correlation 211');
clim([-1,1]);

%% Only take the lower triangle of the connectivity matrix.

triangle111 = tril(con_matrix111, -1);
triangle112 = tril(con_matrix112, -1);
triangle121 = tril(con_matrix121, -1); 
triangle122 = tril(con_matrix122, -1);
triangle211 = tril(con_matrix211, -1);
triangle212 = tril(con_matrix212, -1);

%% Vectorize this lower triangle

triangle111 = triangle111(:);
triangle112 = triangle112(:);
triangle121 = triangle121(:);
triangle122 = triangle122(:);
triangle211 = triangle211(:);
triangle212 = triangle212(:);

%% Compute the AEC. 


% Just a small test, this should be = 1.
aec_value0 = corr2(triangle111, triangle111);

% Just between sessions. Test-restest bestween two different sessions. 
aec_value1 = corr2(triangle111, triangle121);

% Only changing the head to see the influence of the head model.
aec_value2 = corr2(triangle111, triangle112);

% Only changing the sensor data but using the same headmodel. 
aec_value3 = corr2(triangle111, triangle211);

% Changing head model and session. 
aec_value4 = corr2(triangle111, triangle122);


%% Conclusions (if the code is correct of course).



