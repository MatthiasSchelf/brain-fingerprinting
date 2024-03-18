function dasilva = dasilvasensors_PCA(results_PCA, participants)
  
    timepoint1_PCA = {};
    timepoint2_PCA = {};
    
    all_heads = fieldnames(results_PCA);
    for h = 1:length(all_heads)
        current_head = all_heads{h};
        all_timepoints = fieldnames(results_PCA.(current_head));
    
        for t = 1:length(all_timepoints)
            current_timepoint = all_timepoints{t};
            all_forms = fieldnames(results_PCA.(current_head).(current_timepoint));
    
            for f = 1:length(all_forms)
                current_form = all_forms{f};
                end_result = results_PCA.(current_head).(current_timepoint).(current_form);
    
                % Check if the current form is based on sensors
                if contains(lower(current_form), 'sensors') && contains(lower(current_form), '_pca')
                    % Check if the current_timepoint contains 'T1'
                    if contains(lower(current_timepoint), 't1')
                        timepoint1_PCA{end+1} = end_result;
                    % Check if the current_timepoint contains 'T2'
                    elseif contains(lower(current_timepoint), 't2')
                        timepoint2_PCA{end+1} = end_result;
                    end
                end
            end
        end
    end
    
        % Make identifiability matrix
        % Convert cell arrays to matrices for timepoint1 and timepoint2
        timepoint1_matrix_PCA = cell2mat(timepoint1_PCA);
        timepoint2_matrix_PCA = cell2mat(timepoint2_PCA);
        
        % Create an identifiability matrix based on the number of participants
        numParticipants = length(participants);
        Identifiability_matrix_PCA = zeros(numParticipants, numParticipants);
        
        % Calculate the correlation for identifiability matrix
        for i = 1:numParticipants
            for j = 1:numParticipants
                Identifiability_matrix_PCA(i, j) = corr(timepoint1_matrix_PCA(i, :)', timepoint2_matrix_PCA(j, :)');
            end
        end
        
        resultStrings = strings(1, length(participants)); % Initialize a string array to store the results
        
        for p = 1:length(participants)
            % Calculate the differentiability of a participant
            
            % Take the correlation of the participant with itself over time
            selfcorr = Identifiability_matrix_PCA(p, p);
            
            % Take the mean of the correlation of the participant with the others
            column_elements = Identifiability_matrix_PCA(:, p);
            column_elements_without_diagonal = column_elements([1:p-1, p+1:end]);
            meanothercorr = mean(column_elements_without_diagonal);
            
            % Calculate the empirical standard deviation of inter-individual features correlations
            triangle_identifiability_matrix = tril(Identifiability_matrix_PCA, -1);
            
            % Calculate the standard deviation of inter-individual feature correlations
            sd = std(triangle_identifiability_matrix(:));
            
            Dself = (selfcorr - meanothercorr) / sd;
            
        end
        
        % Calculate the mean and median of Dself values
        mean_Dself = mean(Dself);
        median_Dself = median(Dself);
        
        % Display mean and median
        disp(['Mean Dself: ' num2str(mean_Dself)]);
        disp(['Median Dself: ' num2str(median_Dself)]);
        
        % Return the result strings
        dasilva = resultStrings;
end

