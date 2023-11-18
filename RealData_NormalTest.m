
% Real data test for functional connectivity and brain fingerprinting

% Start by opening brainstorm and loading in the necessary data
% Three participants + atlas 

%%

brainstorm

%% Participant 1

%% First check the NaN values

an_indices_Gain = find(isnan(Gain));
disp(['Number of NaN values in Gain: ', num2str(length(Gain))]);

an_indices_GridOrient = find(isnan(GridOrient));
disp(['Number of NaN values in GridOrient: ', num2str(length(GridOrient))]);

%%


%Start by constrianing the data from 3 dimensions into 1

constrained_par1 = bst_gain_orient(Gain, GridOrient);

% Now we have constrained_par1 wich contains 338 sensors x 15002 points (forms the kernel) 


%%
% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 338 senors * 338 sensors x 9001 timepoints

Data_par1_T1 = constrained_par1' * F;

% Now we have 15002 points x 9001 timepoints

%% Altas + Functional connectivity Matrix

%Run over the atlas and remove the NaN

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par1_T1(:,i)=filtfilt(b,a,mean(Data_par1_T1(Scouts.vertices(i),200:npoints),1));
end

% Functional connectivity

con_matrix_par1_T1 = corr(Atlas_par1_T1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1);
colormap("parula");
colorbar;
title('Connectivity Matrix');
clim([-1,1]);






