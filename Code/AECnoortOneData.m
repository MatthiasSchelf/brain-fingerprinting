function AECnoort = AECnoortOneData(participants, fq, fixed_data_T1, fixed_data_T2)
    
    % Initialize a cell array to store results
    resultsCell = cell(length(participants), 1);
    
    % Code for conenctiviy in time domain
    
    parfor p = 1:length(participants)
        
        currentparticipant = participants{p};
        currenthead = strcat(currentparticipant, '_Head');
        % Initialize 'finished' matrix for each participant
        nregions = 68;
        npoints = 9001;
        finished = zeros(npoints - 199, nregions);
        % Store the results in a cell array
        resultsCell{p} = struct('participant', currentparticipant, 'head', currenthead, 'data', struct());
    
        for d = 1:2 %Because there is only timepoint T1 and T2
            
            data_T1 = load(fixed_data_T1)
            data_T2 = load(fixed_data_T2)

            if d == 1
                F = data_T1.F;
            else
                F = data_T2.F;
            end
            
            %decide headmodel, change to the ones you want
            head = currenthead;
    
            %load in the head, data and atlas 
            headData = load(head);
            atlasData = load("scout_Desikan-Killiany_68.mat");
            Gain = headData.Gain;       
            GridOrient = headData.GridOrient;            
        
            % Create a logical index to exclude rows 307 to 339
            rows_to_remove = 307:339;
            
            % Create a logical index for rows to keep for F
            rows_to_keep_F = true(size(F, 1), 1);
            rows_to_keep_F(rows_to_remove) = false;
            F = F(rows_to_keep_F, :);

            % Create a logical index for rows to keep for Gain
            rows_to_keep_Gain = true(size(Gain, 1), 1);
            rows_to_keep_Gain(rows_to_remove) = false;
            Gain = Gain(rows_to_keep_Gain, :);
            
            %Start by constraining the data from 3 dimensions into 1
            
            constrained = bst_gain_orient(Gain, GridOrient);
            
            % Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
            % x timepoints. 
            
            %First transpose so that the time domain is on the first axis.
            F = F';
            
            % Now we filter the data 
            
            fs = 300;
            
            [b, a] = butter(3, fq / (fs/ 2), 'bandpass');  % band to indicate bandpass
            
            filtered = filtfilt(b, a, F);
            
            % Hilbert transform
            
            hil_filtered = abs(hilbert(filtered));
            
            % Now we have  306 sensors x 15002 points (forms the kernel) 
            
            % Combine the constrained data and the timepoints 
            % 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints
            
            Combination = constrained' * hil_filtered' ;
            
            % Now we have 15002 points x 9001 timepoints
            
            % Turn the first field into a cell array. 
            
            % Field to extract
            fieldName = 'Vertices';
            
            % Extract the field and store it in a cell array
            atlas = {atlasData.Scouts.(fieldName)};
            
            % Altas + Functional connectivity Matrix
            
            %Run over the atlas
            
            nregions = 68;
            npoints = 9001;
            
            fs=300;
            
            [b,a]=butter(3,[.5 48]/(fs/2));

            
            for i=1:nregions
            
                finished(:,i)=filtfilt(b,a,mean(Combination(atlas{i},200:npoints),1));
            
            end
            
            % Functional connectivity
            
            con_matrix = corr(finished);
        
            %Indicate what the end result is of this function 
        
            AmplitudeEnvnoort = con_matrix;

            %Now for the sensor data
            currentdata = strcat(currentparticipant, '_T', num2str(d));
            data = currentdata;
            
            % Load in the correct dataset 
            dataData = load(data);
            F = dataData.F; 
            
            % Remove rows 307 until 338 in F this because the recording of these
            % sensors failed + this way all subjects have the same amount of data.
            
            % Create a logical index to exclude rows 307 to 339
            rows_to_remove = 307:339;
            
            % Create a logical index for rows to keep
            rows_to_keep = true(size(F, 1), 1);
            rows_to_keep(rows_to_remove) = false;
            
            % Use logical indexing to remove specified rows
            F = F(rows_to_keep, :);
            
            % Now we filter the data 
            
            %First transpose so that the time domain is on the first axis.
            F = F';
            
            % Now we filter the data 
            
            fs = 300;
            
            [b, a] = butter(3, fq / (fs/ 2), 'bandpass');  % band to indicate bandpass
            
            F_filtered = filtfilt(b, a, F);
            
            % Filter the first 200 timepoints out. 
            
            F_filtered = F_filtered(201:9001, :);
            
            % Hilbert transform
            
            filtered_hil = abs(hilbert(F_filtered));
            
            % Now plot functional connectivity
            
            con_matrix = corr(filtered_hil);
        
            % Indicate what needs to be shown
        
            sensorsAECnoort = con_matrix;
        
            % Store the results in a nested structure
            % Store the results in the cell array
            resultsCell{p}.data.(currentdata).connectivityAEC = AmplitudeEnvnoort;
            resultsCell{p}.data.(currentdata).connectivityAECsensors = sensorsAECnoort;
    
        end
    end

    % Initialize the results structure
    AECnoort = struct();
    
    for p = 1:length(participants)
        current_head = resultsCell{p}.head;
        current_data = fieldnames(resultsCell{p}.data);
    
        % Check if the current_head field exists, if not, create it
        if ~isfield(AECnoort, current_head)
            AECnoort.(current_head) = struct();
        end
    
        % Iterate over the fields in current_data and assign values to the nested structure
        for d = 1:length(current_data)
            current_data_field = current_data{d};
            AECnoort.(current_head).(current_data_field).connectivityAEC = resultsCell{p}.data.(current_data_field).connectivityAEC;
            AECnoort.(current_head).(current_data_field).connectivityAECsensors = resultsCell{p}.data.(current_data_field).connectivityAECsensors;
        end
    end

end