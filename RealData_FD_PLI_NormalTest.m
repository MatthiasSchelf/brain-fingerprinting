
% This code will use the fequency domain and PLI (phase lag index) to
% compute brain fingerprinting

%% first load brainstorm 

brainstorm 

%% load in the atlas

load("scout_Desikan-Killiany_68.mat")

% Participant one at timepoint one

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

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par1_T1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par1_T1 = abs(hilbert(beta_par1_T1));

% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par1_T1 = constrained_par1_T1' * hil_beta_par1_T1';

%% DO PLI for participant 1 at TP 1 

% Parameters
num_points = size(Data_par1_T1, 1);
num_timepoints = size(Data_par1_T1, 2);

% Initialize PLI matrix
PLI_matrix = zeros(num_points, num_points);

% Calculate PLI
for i = 1:num_points
    for j = i+1:num_points
        phase_diff = angle(exp(1i * (Data_par1_T1(j, :) - Data_par1_T1(i, :))));
        PLI_matrix(i, j) = sum(sign(sin(phase_diff))) / num_timepoints;
        PLI_matrix(j, i) = PLI_matrix(i, j); % Since PLI is symmetric
    end
end


%% Now we take the 15002 x 15002 matrix which consists of all the PLI and let it run through the atlas

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

for i=1:nregions

    Atlas_par1_T1(:,i)=filtfilt(b,a,mean(Data_par1_T1(atlas{i},npoints),1));

end

% Functional connectivity

con_matrix_par1_T1 = corr(Atlas_par1_T1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T1');
clim([-1,1]);



















