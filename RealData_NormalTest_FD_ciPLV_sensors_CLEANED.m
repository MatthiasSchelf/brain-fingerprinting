% Here we will use ciPLV to check on the sensors directly. 

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par1_T1 = hilbert(F_filtered);

hil_beta_par1_T1 = angle(hil_beta_par1_T1);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par1_T1');
ndat = hil_beta_par1_T1' ./ abs(hil_beta_par1_T1');
ciplv_par1_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par1_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par1_T1 = corr(ciplv_par1_T1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T1');
clim([-1,1]);

% Participant 1 timepoint 35

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par1_T35 = hilbert(F_filtered);

hil_beta_par1_T35 = angle(hil_beta_par1_T35);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par1_T35');
ndat = hil_beta_par1_T35' ./ abs(hil_beta_par1_T35');
ciplv_par1_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par1_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par1_T35 = corr(ciplv_par1_T35);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T35');
clim([-1,1]);

% Participant 2 timepoint 1

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par2_T1 = hilbert(F_filtered);

hil_beta_par2_T1 = angle(hil_beta_par2_T1);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par2_T1');
ndat = hil_beta_par2_T1' ./ abs(hil_beta_par2_T1');
ciplv_par2_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par2_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par2_T1 = corr(ciplv_par2_T1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T1');
clim([-1,1]);

% Participant 2 timepoint 35

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par2_T35 = hilbert(F_filtered);

hil_beta_par2_T35 = angle(hil_beta_par2_T35);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par2_T35');
ndat = hil_beta_par2_T35' ./ abs(hil_beta_par2_T35');
ciplv_par2_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par2_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par2_T35 = corr(ciplv_par2_T35);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T35');
clim([-1,1]);

% Participant 3 timepoint 1

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par3_T1 = hilbert(F_filtered);

hil_beta_par3_T1 = angle(hil_beta_par3_T1);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par3_T1');
ndat = hil_beta_par3_T1' ./ abs(hil_beta_par3_T1');
ciplv_par3_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par3_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par3_T1 = corr(ciplv_par3_T1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T1');
clim([-1,1]);

% Participant 3 timepoint 35

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

% Now we filter the data so only the beta band remains 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

F_filtered = filtfilt(b, a, F);

% Filter the first 200 timepoints out. 

F_filtered = F_filtered(201:9001, :);

% Hilbert transform

hil_beta_par3_T35 = hilbert(F_filtered);

hil_beta_par3_T35 = angle(hil_beta_par3_T35);

% Calculate ciPLV

[ nc, ns, nt ] = size(hil_beta_par3_T35');
ndat = hil_beta_par3_T35' ./ abs(hil_beta_par3_T35');
ciplv_par3_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par3_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

% Now plot functional connectivity

con_matrix_par3_T35 = corr(ciplv_par3_T35);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T35);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T35');
clim([-1,1]);

%% Now try the PCA  reconstruction 

% Take the lower triangle of all matrices
PCA_triangle_par1_T1 = tril(con_matrix_par1_T1, -1);
PCA_triangle_par1_T35 = tril(con_matrix_par1_T35, -1);
PCA_triangle_par2_T1 = tril(con_matrix_par2_T1, -1);
PCA_triangle_par2_T35 = tril(con_matrix_par2_T35, -1);
PCA_triangle_par3_T1 = tril(con_matrix_par3_T1, -1);
PCA_triangle_par3_T35 = tril(con_matrix_par3_T35, -1);

% Vecrtorize the lower triangle
PCA_triangle_par1_T1 = PCA_triangle_par1_T1(:);
PCA_triangle_par1_T35 = PCA_triangle_par1_T35(:);
PCA_triangle_par2_T1 = PCA_triangle_par2_T1(:);
PCA_triangle_par2_T35 = PCA_triangle_par2_T35(:);
PCA_triangle_par3_T1 = PCA_triangle_par3_T1(:);
PCA_triangle_par3_T35 = PCA_triangle_par3_T35(:);

% Now we want to combine them all. Columns = sessions/participants and rows
% = activity 

PCA_vector_matrix = [PCA_triangle_par1_T1,PCA_triangle_par1_T35,PCA_triangle_par2_T1,PCA_triangle_par2_T35, PCA_triangle_par3_T1,PCA_triangle_par3_T35];

% Now we center the matrix 
mean_PCA_vector_matrix = mean(PCA_vector_matrix);
PCA_vector_matrix = PCA_vector_matrix-mean_PCA_vector_matrix;

% Perform PCA
[coeff, score, latent] = pca(PCA_vector_matrix);

%Select the number of components. 
num_components = 3; % As I do not know the perfect number in this case, I tried them all and this number gave me the highest Idiff
selected_components = coeff(:, 1:num_components);

% Reconstruct the matrix using the selected components
final_PCA_matrix = score(:, 1:num_components) * selected_components';

%% The Amico method 

% Make a cell for the two different timepoints 
timepoint1 = {final_PCA_matrix(:, 1),final_PCA_matrix(:, 3),final_PCA_matrix(:, 5)};
timepoint35 = {final_PCA_matrix(:, 2),final_PCA_matrix(:, 4),final_PCA_matrix(:, 6)};

% Make identifiability matrix
Identifiability_matrix = zeros(3, 3);

for i = 1:3 
    for j = 1:3 
        % Calculate Pearson correlation
        Identifiability_matrix(i, j) = corr2(timepoint1{i}, timepoint35{j});
    end
end

% Visualize the Identifiability matrix
figure;
imagesc(Identifiability_matrix);
clim([0, 1]);
colormap("parula");
colorbar;
title('Identifiability Matrix - Amico Method');
xlabel('Timepoint 35');
ylabel('Timepoint 1');

% Decide Iself Iothers and Idiff

Iself = mean(diag(Identifiability_matrix));
disp(['The Iself of this test example is ', num2str(Iself)]);

triangle_identifiability_matrix = tril(Identifiability_matrix, -1);
triangle_identifiability_matrix = nonzeros(triangle_identifiability_matrix);
Iothers = mean(triangle_identifiability_matrix(:));
disp(['The Iothers of this test example is ', num2str(Iothers)]);

Idiff = (Iself-Iothers)*100;
disp(['The Idiff of this test example is ', num2str(Idiff)]);

%% Da Silva method 

% Make a cell for the two different timepoints 
timepoint1 = {final_PCA_matrix(:, 1),final_PCA_matrix(:, 3),final_PCA_matrix(:, 5)};
timepoint35 = {final_PCA_matrix(:, 2),final_PCA_matrix(:, 4),final_PCA_matrix(:, 6)};

% Make identifiability matrix
Identifiability_matrix = zeros(3, 3);

for i = 1:3 
    for j = 1:3 
        % Calculate Pearson correlation
        Identifiability_matrix(i, j) = corr2(timepoint1{i}, timepoint35{j});
    end
end

% Calculate the differentiability of participant one 

% Take the correlation of participant one with itself over time
selfcorrpar1 = Identifiability_matrix(1,1);

% Take the mean of the correlation of participant one with the others
column_index_par1 = 1;
column_elements_par1 = Identifiability_matrix(:, column_index_par1);
column_elements_without_diagonal_par1 = column_elements_par1([1:column_index_par1-1, column_index_par1+1:end]);
meanothercorrpar1 = mean(column_elements_without_diagonal_par1);

% calculate the empirical standard deviation of inter-individual features correlations
triangle_identifiability_matrix_par1 = tril(Identifiability_matrix, -1);
% Calculate the standard deviation of inter-individual feature correlations
sdpar1 = std(triangle_identifiability_matrix_par1(:));

Dselfpar1 = (selfcorrpar1 - meanothercorrpar1)/sdpar1;
disp(['The Dself of participant one is ', num2str(Dselfpar1)]);

% Calculate the differentiability of participant two 

% Take the correlation of participant one with itself over time
selfcorrpar2 = Identifiability_matrix(2,2);

% Take the mean of the correlation of participant one with the others
column_index_par2 = 2;
column_elements_par2 = Identifiability_matrix(:, column_index_par2);
column_elements_without_diagonal_par2 = column_elements_par2([1:column_index_par2-1, column_index_par2+1:end]);
meanothercorrpar2 = mean(column_elements_without_diagonal_par2);

% calculate the empirical standard deviation of inter-individual features correlations
triangle_identifiability_matrix_par2 = tril(Identifiability_matrix, -1);
% Calculate the standard deviation of inter-individual feature correlations
sdpar2 = std(triangle_identifiability_matrix_par2(:));

Dselfpar2 = (selfcorrpar2 - meanothercorrpar2)/sdpar2;
disp(['The Dself of participant two is ', num2str(Dselfpar2)]);

% Calculate the differentiability of participant three

% Take the correlation of participant one with itself over time
selfcorrpar3 = Identifiability_matrix(3,3);

% Take the mean of the correlation of participant one with the others
column_index_par3 = 3;
column_elements_par3 = Identifiability_matrix(:, column_index_par3);
column_elements_without_diagonal_par3 = column_elements_par3([1:column_index_par3-1, column_index_par3+1:end]);
meanothercorrpar3 = mean(column_elements_without_diagonal_par3);

% calculate the empirical standard deviation of inter-individual features correlations
triangle_identifiability_matrix_par3 = tril(Identifiability_matrix, -1);
% Calculate the standard deviation of inter-individual feature correlations
sdpar3 = std(triangle_identifiability_matrix_par3(:));

Dselfpar3 = (selfcorrpar3 - meanothercorrpar3)/sdpar3;
disp(['The Dself of participant three is ', num2str(Dselfpar3)]);


