
function AECort = AECortMatched(participants, fq)
    
    % Initialize a cell array to store results
    resultsCell = cell(length(participants), 1);
        
        % Code for conenctiviy in time domain
        
        parfor p = 1:length(participants)
            
            currentparticipant = participants{p};
            currenthead = strcat(currentparticipant, '_Head');
            % Store the results in a cell array
            resultsCell{p} = struct('participant', currentparticipant, 'head', currenthead, 'data', struct());
        
            for d = 1:2 %Because there is only timepoint T1 and T2
                
                currentdata = strcat(currentparticipant, '_T', num2str(d));
                
                %decide headmodel, change to the ones you want
                head = currenthead;
            
                %decide data 
                data = currentdata;
            
                %Load in the necessary data 
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
                
                % Now we filter the data 
                
                fs = 300;
                
                [b, a] = butter(3, fq / (fs/ 2), 'bandpass');  % band to indicate bandpass
                
                filtered = filtfilt(b, a, F);
                
                % AEC with orthogonalization by using function in script RealData_AEC_Orthogonalization
                
                % Now we call the function and let it work 
                
                a = filtered;
                
                filtered_AEC = AEC(a);
                
                % Get the full dataset
                
                Data = constrained' *filtered_AEC' ;
                
                % Go through the atlas 
                
                % Field to extract
                fieldName = 'Vertices';
                
                % Extract the field and store it in a cell array
                atlas = {atlasData.Scouts.(fieldName)};
                
                % Altas + Functional connectivity Matrix
                
                %Run over the atlas
                
                nregions = 68;
                npoints = 306;
                
                fs=300;
                
                [b,a]=butter(3,[.5 48]/(fs/2));
                
                % Initialize Atlas_par1_T1 matrix
                finished = zeros(npoints, nregions);
                
                for i = 1:nregions
                    region_data = Data(atlas{i}, 1:npoints);
                
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
                
                % Functional connectivity
                
                con_matrix = corr(finished);
            
                %Indicate what the result needs to be 
            
                connectivityAECort = con_matrix;

                % Load in correct dataset 
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
                
                % AEC with orthogonalization by using function in script RealData_AEC_Orthogonalization
                
                % Now we call the function and let it work 
                
                a = F_filtered;
                
                AECorth = AEC(a);
                
                % Now plot functional connectivity
                
                con_matrix = corr(AECorth);
            
                %Indicate what the result is 
            
                connectivityAECortsensors = con_matrix;
            
                % Store the results in a nested structure
                resultsCell{p}.data.(currentdata).connectivityAECort = connectivityAECort;
                resultsCell{p}.data.(currentdata).connectivityAECortsensors = connectivityAECortsensors;
        
            end
        end

        % Initialize the results structure
        AECort = struct();
        
        for p = 1:length(participants)
            current_head = resultsCell{p}.head;
            current_data = fieldnames(resultsCell{p}.data);
        
            % Check if the current_head field exists, if not, create it
            if ~isfield(AECort, current_head)
                AECort.(current_head) = struct();
            end
        
            % Iterate over the fields in current_data and assign values to the nested structure
            for d = 1:length(current_data)
                current_data_field = current_data{d};
                AECort.(current_head).(current_data_field).connectivityAECort = resultsCell{p}.data.(current_data_field).connectivityAECort;
                AECort.(current_head).(current_data_field).connectivityAECortsensors = resultsCell{p}.data.(current_data_field).connectivityAECortsensors;
            end
        end

end