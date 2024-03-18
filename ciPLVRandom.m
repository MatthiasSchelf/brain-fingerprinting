function ciPLV =ciPLVRandom(participants, fq, shuffledheads)
    
    % Shuffle the heads outside the parfor loop
    shuffledOrder = randperm(length(participants));
    shuffledheads = shuffledheads(shuffledOrder);
    
    % Initialize a cell array to store results
    resultsCell = cell(length(participants), 1);
    
    % Code for conenctiviy in time domain
    
    for p = 1:length(participants)
        
        currentparticipant = participants{p};
        currenthead = strcat(shuffledheads{p}, '_Head');
    
        for d = 1:2 %Because there is only timepoint T1 and T2
            
            currentdata = strcat(currentparticipant, '_T', num2str(d));
            
            %decide headmodel, change to the ones you want
            head = currenthead;
        
            %decide data 
            data = currentdata;
        
             %Load in the head, data and atlas. 
        
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
            
            % Now we need to filter to the beta frequency in the F matrix. As the F matrix has the sensors
            % x timepoints. 
            
            %First transpose so that the time domain is on the first axis.
            F = F';
            
            % Now we filter the data so only the beta waves remain
            
            fs = 300;
            
            [b, a] = butter(3, fq / (fs/ 2), 'bandpass');  % band to indicate bandpass
            
            filtered = filtfilt(b, a, F);
            
            % Hilbert transform
            
            hil_filtered = hilbert(filtered);
            
            hil_angle_filtered = angle(hil_filtered);
            
            % Get the full dataset
            
            Data = constrained' *hil_angle_filtered' ;
            
            % Calculate ciPLV
            
            [ nc, ns, nt ] = size(Data);
            ndat = Data ./ abs(Data);
            ciplv = zeros(nc, nc, nt);
            
            for t = 1:nt
                cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
                
                % Calculate PLV
                plv = abs(cross_spectrum) / ns;
                
                % Calculate ciPLV by subtracting the imaginary part
                ciplv(:, :, t) = plv - imag(cross_spectrum) / ns;
            end
            
            %
            % Now we need to run over the atlas. 
            
            % Now we have 15002 points x 15002 ciPLV's
            
            % Turn the first field into a cell array. 
            
            % Field to extract
            fieldName = 'Vertices';
            
            % Extract the field and store it in a cell array
            atlas = {atlasData.Scouts.(fieldName)};
            
            % Altas + Functional connectivity Matrix
            
            %Run over the atlas
            
            nregions = 68;
            npoints = 15002;
            
            fs=300;
            
            [b,a]=butter(3,[.5 48]/(fs/2));
            
            % Initialize Atlas_par1_T1 matrix
            finished = zeros(npoints, nregions);
            
            for i = 1:nregions
                region_data = ciplv(atlas{i}, 1:npoints);
            
                % Check if the length of data is sufficient for filtering
                if size(region_data, 1) >= 18
                    filtered_data = filtfilt(b, a, mean(region_data, 1));
                    
                    % Ensure that the size of filtered_data matches the size of Atlas_par1_T1
                    finished(1:numel(filtered_data), i) = filtered_data(:);
                else
                    % Handle the case when the data length is insufficient
                    disp(['Skipping region ', num2str(i), ' due to insufficient data length.']);
                end
            end
            
            %
            
            % Functional connectivity
            
            con_matrix = corr(finished);
            
            %Indicate what the end result is of this function 
        
            connectivityciPLV = con_matrix;
        
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
            
            % Now we filter the data so only the beta band remains 
            
            fs = 300;
            
            [b, a] = butter(3, fq / (fs/ 2), 'bandpass');  % band to indicate bandpass
            
            F_filtered = filtfilt(b, a, F);
            
            % Filter the first 200 timepoints out. 
            
            F_filtered = F_filtered(201:9001, :);
            
            % Hilbert transform
            
            hil_filtered = hilbert(F_filtered);
            
            hil_filtered = angle(hil_filtered);
            
            % Calculate ciPLV
            
            [ nc, ns, nt ] = size(hil_filtered');
            ndat = hil_filtered' ./ abs(hil_filtered');
            ciplv = zeros(nc, nc, nt);
            
            for t = 1:nt
                cross_spectrum = ndat(:, :, t) * ndat(:, :, t)';
                
                % Calculate PLV
                plv = abs(cross_spectrum) / ns;
                
                % Calculate ciPLV by subtracting the imaginary part
                ciplv(:, :, t) = plv - imag(cross_spectrum) / ns;
            end
            
            % Now plot functional connectivity
            
            con_matrix = corr(ciplv);
        
            %Indicate what the result is 
        
            sensorciPLV = con_matrix;

            % Store the head information in resultsCell
            resultsCell{p}.head = currenthead;
        
            % Store the results in a nested structure
            resultsCell{p}.data.(currentdata).connectivityciPLV = connectivityciPLV;
            resultsCell{p}.data.(currentdata).connectivityciPLVsensors = sensorciPLV;
    
        end
    end
    
    % Initialize the results structure
    ciPLV = struct();
    
    for p = 1:length(participants)
        current_head = resultsCell{p}.head;
        current_data = fieldnames(resultsCell{p}.data);
    
        % Check if the current_head field exists, if not, create it
        if ~isfield(ciPLV, current_head)
            ciPLV.(current_head) = struct();
        end
    
        % Iterate over the fields in current_data and assign values to the nested structure
        for d = 1:length(current_data)
            current_data_field = current_data{d};
            ciPLV.(current_head).(current_data_field).connectivityciPLV = resultsCell{p}.data.(current_data_field).connectivityciPLV;
            ciPLV.(current_head).(current_data_field).connectivityciPLVsensors = resultsCell{p}.data.(current_data_field).connectivityciPLVsensors;
        end
    end

end