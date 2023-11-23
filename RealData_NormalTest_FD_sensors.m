
% Within this dataset fingerprinting is done with just using the sensors
% and this within the beta band. 

%% Participant 1 timepoint 1 

% Load in correct dataset 

load("Par1_T1.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par1_T1 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T1');
clim([-1,1]);

%% Participant 1 timepoint 35

% Load in correct dataset 

load("Par1_T35.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par1_T35 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T35');
clim([-1,1]);

%% Participant 2 timepoint 1

% Load in correct dataset 

load("Par2_T1.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par2_T1 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T1');
clim([-1,1]);

%% Participant 2 timepoint 35

% Load in correct dataset 

load("Par2_T35.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par2_T35 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T35');
clim([-1,1]);

%% Participant 3 timepoint 1

% Load in correct dataset 

load("Par3_T1.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par3_T1 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T1');
clim([-1,1]);

%% Participant 3 timepoint 35

% Load in correct dataset 

load("Par3_T35.mat");

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Now we filter the data 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_F = abs(hilbert(F_filtered));

% Now plot functional connectivity

con_matrix_par3_T35 = corr(hil_F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T35');
clim([-1,1]);


%% Now let's look at the corerlations

% Take the lower triangle under the diagonal.

triangle_par1_T1 = tril(con_matrix_par1_T1, -1);
triangle_par1_T35 = tril(con_matrix_par1_T35, -1);
triangle_par2_T1 = tril(con_matrix_par2_T1, -1); 
triangle_par2_T35 = tril(con_matrix_par2_T35, -1);
triangle_par3_T1 = tril(con_matrix_par3_T1, -1);
triangle_par3_T35 = tril(con_matrix_par3_T35, -1);


% Now we vectorize these triangles

triangle_par1_T1 = triangle_par1_T1(:);
triangle_par1_T35 = triangle_par1_T35(:);
triangle_par2_T1 = triangle_par2_T1(:);
triangle_par2_T35 = triangle_par2_T35(:);
triangle_par3_T1 = triangle_par3_T1(:);
triangle_par3_T35 = triangle_par3_T35(:);

% Correlations calculations

corr_par1_T1_T35 = corr(triangle_par1_T1, triangle_par1_T35);
disp(['The correlation between T1 and T35 of participant one is ', num2str(corr_par1_T1_T35)]);

corr_par2_T1_T35 = corr(triangle_par2_T1, triangle_par2_T35);
disp(['The correlation between T1 and T35 of participant two is ', num2str(corr_par2_T1_T35)]);

corr_par2_T1_T35 = corr(triangle_par3_T1, triangle_par3_T35);
disp(['The correlation between T1 and T35 of participant three is ',num2str(corr_par2_T1_T35)]);

% Other calculations

corr_par1_par2 = corr(triangle_par1_T1, triangle_par2_T35);
disp(['The correlation between par1 and par2 (from T1 to T35) is ', num2str(corr_par1_par2)]);

corr_par1_par3 = corr(triangle_par1_T1, triangle_par3_T35);
disp(['The correlation between par1 and par3 (from T1 to T35) is ', num2str(corr_par1_par3)]);

corr_par2_par1 = corr(triangle_par2_T1, triangle_par1_T35);
disp(['The correlation between par2 and par1 (from T1 to T35) is ', num2str(corr_par2_par1)]);

corr_par2_par3 = corr(triangle_par2_T1, triangle_par3_T35);
disp(['The correlation between par2 and par3 (from T1 to T35) is ', num2str(corr_par2_par3)]);

corr_par3_par1 = corr(triangle_par3_T1, triangle_par1_T35);
disp(['The correlation between par3 and par1 (from T1 to T35) is ', num2str(corr_par3_par1)]);

corr_par3_par2 = corr(triangle_par3_T1, triangle_par2_T35);
disp(['The correlation between par3 and par2 (from T1 to T35) is ', num2str(corr_par3_par2)]);

