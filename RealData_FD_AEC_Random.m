
% Witin this dataset, the Frequency Domain will be investigated and random
% head models will be used. To see which head models are used when, check
% the PDF document. 

%%

brainstorm

%%
load('scout_Desikan-Killiany_68.mat')

%% Participant 1 at timepoint 1 with head 3

%load data
load("headmodel_Par3.mat")
load("Par1_T1.mat")

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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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



% Participant 1 at timepoint 35 and head 2 

%load in data
load("headmodel_Par2.mat")
load("Par1_T35.mat")


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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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


% Participant 2 at timepoint 1 and head 1
load("headmodel_Par1.mat")
load("Par2_T1.mat")


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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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


% Participant 2 at timepoint 35 with head 3
load("headmodel_Par3.mat")
load("Par2_T35.mat")

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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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

% Participant 3 at timepoint 1 and head 2
load("headmodel_Par2.mat")
load("Par3_T1.mat")

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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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


% Participant 3 at timepoint 35 and head 1

%Load in data 
load("headmodel_Par1.mat")
load("Par3_T35.mat")

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

[b, a] = butter(3, [0.5,4] / (fs/ 2), 'bandpass');  % band to indicate bandpass

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


%Now calculate the correlations. 

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

%% Amico method 

% Make a cell for the two different timepoints 
timepoint1 = {triangle_par1_T1,triangle_par2_T1,triangle_par3_T1};
timepoint35 = {triangle_par1_T35,triangle_par2_T35,triangle_par3_T35};

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
timepoint1 = {triangle_par1_T1,triangle_par2_T1,triangle_par3_T1};
timepoint35 = {triangle_par1_T35,triangle_par2_T35,triangle_par3_T35};

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











