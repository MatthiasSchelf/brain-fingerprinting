
%Transpose first so that the matrix is time X channels (as indicated in the
%AEC script)
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

% First parameters
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


%% Once filtered, use the AEC function
a = beta_sub01ses01;
N=size(a,2);
AEC(1:N,1:N)=0;
complex_a=hilbert(a);

for i=1:N
    for j=1:N
        if i<j        
        ort1=orthog_timedomain(a(:,i),a(:,j));
        complex_ort1=abs(hilbert(ort1));
        AEC1=abs(corrcoef(complex_ort1,abs(complex_a(:,i))));        
        ort2=orthog_timedomain(a(:,j),a(:,i));
        complex_ort2=abs(hilbert(ort2));
        AEC2=abs(corrcoef(complex_ort2,abs(complex_a(:,j))));        
        AEC_mean=(AEC1(1,2)+AEC2(1,2))/2;        
        AEC(i,j)=AEC_mean;       
        AEC(j,i)=AEC(i,j);        
        end        
    end
end
