
% Witin this dataset, the Frequency Domain will be investigated and random
% head models will be used. To see which head models are used when, check
% the PDF document. 

%%

brainstorm

%% Participant 1 at timepoint 1 with head 3

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par1_T1_head3 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par1_T1_head3 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par1_T1_head3 = abs(hilbert(beta_par1_T1_head3));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par1_T1_head3 = constrained_par1_T1_head3' * hil_beta_par1_T1_head3' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par1_T1_head3(:,i)=filtfilt(b,a,mean(Data_par1_T1_head3(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par1_T1_head3 = corr(Atlas_par1_T1_head3);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T1_head3);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T1 head3');
clim([-1,1]);

%% Participant 1 at timepoint 35 and head 2 

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par1_T35_head2 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par1_T35_head2 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par1_T35_head2 = abs(hilbert(beta_par1_T35_head2));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par1_T35_head2 = constrained_par1_T35_head2' * hil_beta_par1_T35_head2' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par1_T35_head2(:,i)=filtfilt(b,a,mean(Data_par1_T35_head2(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par1_T35_head2 = corr(Atlas_par1_T35_head2);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par1_T35_head2);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 1 T35 head2');
clim([-1,1]);

%% Participant 2 at timepoint 1 and head 1

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par2_T1_head1 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par2_T1_head1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par2_T1_head1 = abs(hilbert(beta_par2_T1_head1));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par2_T1_head1 = constrained_par2_T1_head1' * hil_beta_par2_T1_head1' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par2_T1_head1(:,i)=filtfilt(b,a,mean(Data_par2_T1_head1(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par2_T1_head1 = corr(Atlas_par2_T1_head1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T1_head1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T1 head 1');
clim([-1,1]);

%% Participant 2 at timepoint 35 with head 3

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par2_T35_head3 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par2_T35_head3 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par2_T35_head3 = abs(hilbert(beta_par2_T35_head3));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par2_T35_head3 = constrained_par2_T35_head3' * hil_beta_par2_T35_head3' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par2_T35_head3(:,i)=filtfilt(b,a,mean(Data_par2_T35_head3(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par2_T35_head3 = corr(Atlas_par2_T35_head3);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par2_T35_head3);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 2 T35 head3');
clim([-1,1]);

%% Participant 3 at timepoint 1 and head 2

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par3_T1_head2 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par3_T1_head2 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par3_T1_head2 = abs(hilbert(beta_par3_T1_head2));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par3_T1_head2 = constrained_par3_T1_head2' * hil_beta_par3_T1_head2' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par3_T1_head2(:,i)=filtfilt(b,a,mean(Data_par3_T1_head2(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par3_T1_head2 = corr(Atlas_par3_T1_head2);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T1_head2);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T1 head 2');
clim([-1,1]);

%% Participant 3 at timepoint 35 and head 1

% Remove rows 307 until 338 in Gain and F

% Define row indices to remove
indices_to_remove = [307:338];

% Remove rows based on indices in Gain and F
Gain(indices_to_remove, :) = [];

F(indices_to_remove, :) = [];

%Start by constraining the data from 3 dimensions into 1

constrained_par3_T35_head1 = bst_gain_orient(Gain, GridOrient);

% Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
% x timepoints. 

%First transpose so that the time domain is on the first axis.
F = F';

% Now we filter the data 

fs = 300;

[b, a] = butter(3, [13, 30] / (fs/ 2), 'bandpass');  % band to indicate bandpass

beta_par3_T35_head1 = filtfilt(b, a, F);

% Hilbert transform

hil_beta_par3_T35_head1 = abs(hilbert(beta_par3_T35_head1));

% Now we have  306 sensors x 15002 points (forms the kernel) 

% Combine the constrained data and the timepoints 
% 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints

Data_par3_T35_head1 = constrained_par3_T35_head1' * hil_beta_par3_T35_head1' ;

% Now we have 15002 points x 9001 timepoints

% Turn the first field into a cell array. 

% Field to extract
fieldName = 'Vertices';

% Extract the field and store it in a cell array
atlas = {Scouts.(fieldName)};

% Altas + Functional connectivity Matrix

%Run over the atlas

nregions = 68;
npoints = 9001;

fs=300;

[b,a]=butter(3,[.5 48]/(fs/2));

for i=1:nregions

    Atlas_par3_T35_head1(:,i)=filtfilt(b,a,mean(Data_par3_T35_head1(atlas{i},200:npoints),1));

end

% Functional connectivity

con_matrix_par3_T35_head1 = corr(Atlas_par3_T35_head1);

% Plot connectivity matrix

figure;
imagesc(con_matrix_par3_T35_head1);
colormap("parula");
colorbar;
title('Connectivity Matrix Participant 3 T35 head 1');
clim([-1,1]);

%% Now calculate the correlations. 

% Take the lower triangle under the diagonal.

triangle_par1_T1 = tril(con_matrix_par1_T1_head3, -1);
triangle_par1_T35 = tril(con_matrix_par1_T35_head2, -1);
triangle_par2_T1 = tril(con_matrix_par2_T1_head1, -1); 
triangle_par2_T35 = tril(con_matrix_par2_T35_head3, -1);
triangle_par3_T1 = tril(con_matrix_par3_T1_head2, -1);
triangle_par3_T35 = tril(con_matrix_par3_T35_head1, -1);


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

% Other correlations

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







