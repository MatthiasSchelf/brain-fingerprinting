
% Here we test brain fingerprinting solely on the sensors. 

%% Participant 1 at timepoint 1

% Load in the correct dataset 
load("Par1_T1.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par1_T1 = corr(F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T1');
clim([-1,1]);

%% Particiant 1 at timepoint 35 

% Load in the correct dataset 
load("Par1_T35.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

% Filter some 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par1_T35 = corr(F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T35');
clim([-1,1]);

%% Particiant 2 at timepoint 1

% Load in the correct dataset 
load("Par2_T1.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

% Filter some 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par2_T1 = corr(F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T1');
clim([-1,1]);

%% Particiant 2 at timepoint 35

% Load in the correct dataset 
load("Par2_T35.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

% Filter some 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par2_T35 = corr(F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T35');
clim([-1,1]);

%% Particiant 3 at timepoint 1

% Load in the correct dataset 
load("Par3_T1.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

% Filter some 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par3_T1 = corr(F);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T1');
clim([-1,1]);

%% Particiant 3 at timepoint 35

% Load in the correct dataset 
load("Par3_T35.mat")

% Remove rows 307 until 338 in F this because the recording of these
% sensors failed + this way all subjects have the same amount of data.

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in F

F(indices_to_remove, :) = [];

% Filter the data just a little bit 

% Filter some 

npoints = 9001;
fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

F = filtfilt(b, a, F');

% Filter the first 200 timepoints out. 

F_filtered = F(201:9001, :);

% Now that the data are appropiately filtered, construct the connectivity matrix

con_matrix_par3_T35 = corr(F);

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




