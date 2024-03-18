function dasilva = dasilva(results, participants)
  
    timepoint1 = {};
    timepoint2 = {};
    
    all_heads = fieldnames(results);
    for h = 1:length(all_heads)
        current_head = all_heads{h};
        all_timepoints = fieldnames(results.(current_head));
    
        for t = 1:length(all_timepoints)
            current_timepoint = all_timepoints{t};
            all_forms = fieldnames(results.(current_head).(current_timepoint));
    
            for f = 1:length(all_forms)
                current_form = all_forms{f};
                end_result = results.(current_head).(current_timepoint).(current_form);
    
                % Check if the current form is based on sensors
                if ~contains(lower(current_form), 'sensors') && ~contains(lower(current_form), '_pca')
                    % Check if the current_timepoint contains 'T1'
                    if contains(lower(current_timepoint), 't1')
                        timepoint1{end+1} = end_result;
                    % Check if the current_timepoint contains 'T2'
                    elseif contains(lower(current_timepoint), 't2')
                        timepoint2{end+1} = end_result;
                    end
                end
            end
        end
    end
    
    % Make identifiability matrix
    Identifiability_matrix = zeros(length(participants), length(participants));
    
    for i = 1:length(participants) 
        for j = 1:length(participants) 
            % Calculate Pearson correlation
            Identifiability_matrix(i, j) = corr2(timepoint1{i}, timepoint2{j});
        end
    end
    
    resultStrings = strings(1, length(participants)); % Initialize a string array to store the results
    Dself_values = zeros(1, length(participants)); % Initialize an array to store Dself values
    
    for p = 1:length(participants)
        % Calculate the differentiability of a participant
        
        % Take the correlation of the participant with itself over time
        selfcorr = Identifiability_matrix(p, p);
        
        % Take the mean of the correlation of the participant with the others
        column_elements = Identifiability_matrix(:, p);
        column_elements_without_diagonal = column_elements([1:p-1, p+1:end]);
        meanothercorr = mean(column_elements_without_diagonal);
        
        % Calculate the empirical standard deviation of inter-individual features correlations
        triangle_identifiability_matrix = tril(Identifiability_matrix, -1);
        
        % Calculate the standard deviation of inter-individual feature correlations
        sd = std(triangle_identifiability_matrix(:));
        
        Dself = (selfcorr - meanothercorr) / sd;
        
        % Store the result in the arrays
        resultStrings(p) = sprintf('Participant %d: %f', p, Dself);
        Dself_values(p) = Dself;
    end
    
    % Calculate the mean and median of Dself values
    mean_Dself = mean(Dself_values);
    median_Dself = median(Dself_values);
    
    % Display mean and median
    disp(['Mean Dself: ' num2str(mean_Dself)]);
    disp(['Median Dself: ' num2str(median_Dself)]);
    
    % Return the result strings
    dasilva = resultStrings;

end

