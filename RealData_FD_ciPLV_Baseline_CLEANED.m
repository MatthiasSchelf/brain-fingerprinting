%% first load brainstorm 

brainstorm 

%% Load in the atlas

load("scout_Desikan-Killiany_68.mat")

%% Participant one at timepoint one

%load in the data

load("headmodel_Par1.mat")
load("Par1_T1.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par1_T1 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par1_T1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par1_T1 = hilbert(beta_par1_T1);

hil_beta_par1_T1 = angle(hil_beta_par1_T1);

% Get the full dataset

Data_par1_T1 = constrained_par1_T1' *hil_beta_par1_T1' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par1_T1);
ndat = Data_par1_T1 ./ abs(Data_par1_T1);
ciplv_par1_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par1_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par1_T1 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par1_T1(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par1_T1(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

%

% Functional connectivity

con_matrix_par1_T1 = corr(Atlas_par1_T1);


% Participant one at timepoint 35

%load in the data

load("headmodel_Par1.mat")
load("Par1_T35.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par1_T35 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par1_T35 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par1_T35 = hilbert(beta_par1_T35);

hil_beta_par1_T35 = angle(hil_beta_par1_T35);

% Get the full dataset

Data_par1_T35 = constrained_par1_T35' *hil_beta_par1_T35' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par1_T35);
ndat = Data_par1_T35 ./ abs(Data_par1_T35);
ciplv_par1_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par1_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par1_T35 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par1_T35(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par1_T35(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

%

% Functional connectivity

con_matrix_par1_T35 = corr(Atlas_par1_T35);

% Participant 2 at timepoint 1

%load in the data

load("headmodel_Par1.mat")
load("Par2_T1.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par2_T1 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par2_T1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par2_T1 = hilbert(beta_par2_T1);

hil_beta_par2_T1 = angle(hil_beta_par2_T1);

% Get the full dataset

Data_par2_T1 = constrained_par2_T1' *hil_beta_par2_T1' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par1_T35);
ndat = Data_par2_T1 ./ abs(Data_par2_T1);
ciplv_par2_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par2_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par2_T1 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par2_T1(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par2_T1(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

%

% Functional connectivity

con_matrix_par2_T1 = corr(Atlas_par2_T1);

% Participant 2 at timepoint 35

%load in the data

load("headmodel_Par1.mat")
load("Par2_T35.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par2_T35 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par2_T35 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par2_T35 = hilbert(beta_par2_T35);

hil_beta_par2_T35 = angle(hil_beta_par2_T35);

% Get the full dataset

Data_par2_T35 = constrained_par2_T35' *hil_beta_par2_T35' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par2_T35);
ndat = Data_par2_T35 ./ abs(Data_par2_T35);
ciplv_par2_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par2_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par2_T35 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par2_T35(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par2_T35(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

%

% Functional connectivity

con_matrix_par2_T35 = corr(Atlas_par2_T35);

% Participant 3 at timepoint 1

%load in the data

load("headmodel_Par1.mat")
load("Par3_T1.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par3_T1 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par3_T1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par3_T1 = hilbert(beta_par3_T1);

hil_beta_par3_T1 = angle(hil_beta_par3_T1);

% Get the full dataset

Data_par3_T1 = constrained_par3_T1' *hil_beta_par3_T1' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par3_T1);
ndat = Data_par3_T1 ./ abs(Data_par3_T1);
ciplv_par3_T1 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par3_T1(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par3_T1 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par3_T1(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par3_T1(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

%

% Functional connectivity

con_matrix_par3_T1 = corr(Atlas_par3_T1);


% Participant 3 at timepoint 35

%load in the data

load("headmodel_Par1.mat")
load("Par3_T35.mat")

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par3_T35 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data so only the beta waves remain

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par3_T35 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par3_T35 = hilbert(beta_par3_T35);

hil_beta_par3_T35 = angle(hil_beta_par3_T35);

% Get the full dataset

Data_par3_T35 = constrained_par3_T35' *hil_beta_par3_T35' ;

% Calculate ciPLV

[ nc, ns, nt ] = size(Data_par3_T35);
ndat = Data_par3_T35 ./ abs(Data_par3_T35);
ciplv_par3_T35 = zeros(nc, nc, nt);

for t = 1:nt
    cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
    
    % Calculate PLV
    plv = abs(cross_spectrum) / ns;
    
    % Calculate ciPLV by subtracting the imaginary part
    ciplv_par3_T35(:, :, t) = plv - imag(cross_spectrum) / ns;
end

%
% Now we need to run over the atlas. 

% Now we have 15002 points x 15002 ciPLV's

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 15002;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

% Initialize Atlas_par1_T1 matrix
Atlas_par3_T35 = zeros(npoints, nregions);

for i = 1:nregions
    region_data = ciplv_par3_T35(atlas{i}, 1:npoints);

    % Check if the length of data is sufficient for filtering
    if size(region_data, 1) >= 18
        filtered_data = filtfilt(b, a, mean(region_data, 1));
        
        % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
        Atlas_par3_T35(1:numel(filtered_data), i) = filtered_data(:);
    else
        % Handle the case when the data length is insufficient
        disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
    end
end

% Functional connectivity

con_matrix_par3_T35 = corr(Atlas_par3_T35);

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



