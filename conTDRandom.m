function conTD = conTDRandom(participants,shuffledheads)
    
    % Shuffle the heads outside the parfor loop
    shuffledOrder = randperm(length(participants));
    shuffledheads = shuffledheads(shuffledOrder);
    
    % Initialize a cell array to store results
    resultsCell = cell(length(participants), 1);
    
    % Code for conenctiviy in time domain
    
    parfor p = 1:length(participants)
        
        currentparticipant = participants{p};
        currenthead = strcat(shuffledheads{p}, '_Head');
        % Initialize 'finished' matrix for each participant
        nregions = 68;
        npoints = 9001;
        finished = zeros(npoints - 199, nregions);
        % Store the results in a cell array
        resultsCell{p} = struct('participant', currentparticipant, 'head', currenthead, 'data', struct());
    
        for d = 1:2 %Because there is only timepoint T1 and T2
            
            currentdata = strcat(currentparticipant, '_T', num2str(d));
            
            %decide headmodel, change to the ones you want
            head = currenthead;
        
            %decide data 
            data = currentdata;
    
            %load in the head, data and atlas 

            headData = load(head);
            dataData = load(data);
            atlasData = load("scout_Desikan-Killiany_68.mat");
            Gain = headData.Gain;       
            GridOrient = headData.GridOrient; 
            F = dataData.F;    
        
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
            
            % Now we have constrained_par1 wich contains 306 sensors x 15002 points (forms the kernel) 
            
            % Combine the constrained data and the timepoints 
            % 15002 points (forms the kernel) x 306 senors * 306 sensors x 9001 timepoints
            
            Combination = constrained' * F;
            
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
        
            connectivityTD = con_matrix;
        
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
            
            % Filter the data just a little bit 
            
            npoints = 9001;
            fs=300;
            
            [b,a]=butter(3,[.5 48]/(fs/2));
            
            F = filtfilt(b, a, F');
            
            % Filter the first 200 timepoints out. 
            
            F_filtered = F(201:9001, :);
            
            % Now that the data are appropiately filtered, construct the connectivity matrix
            
            con_matrix = corr(F_filtered);
        
            % Indicate what needs to be given back to the source code
            sensorTD = con_matrix;
            
            % Store the results in the cell array
            resultsCell{p}.data.(currentdata).connectivityTD = connectivityTD;
            resultsCell{p}.data.(currentdata).connectivityTDsensors = sensorTD;
    
        end
    end

    % Initialize the results structure
    conTD = struct();
    
    for p = 1:length(participants)
        current_head = resultsCell{p}.head;
        current_data = fieldnames(resultsCell{p}.data);
    
        % Check if the current_head field exists, if not, create it
        if ~isfield(conTD, current_head)
            conTD.(current_head) = struct();
        end
    
        % Iterate over the fields in current_data and assign values to the nested structure
        for d = 1:length(current_data)
            current_data_field = current_data{d};
            conTD.(current_head).(current_data_field).connectivityTD = resultsCell{p}.data.(current_data_field).connectivityTD;
            conTD.(current_head).(current_data_field).connectivityTDsensors = resultsCell{p}.data.(current_data_field).connectivityTDsensors;
        end
    end

end

