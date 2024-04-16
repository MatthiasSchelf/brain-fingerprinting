function dasilva = dasilvasensors(results, participants)
  
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
                if contains(lower(current_form), 'sensors') && ~contains(lower(current_form), '_pca')
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
    
    % Make identifiability matrices
    Pearson_Identifiability_matrix = zeros(length(participants), length(participants));
    covstatis_Identifiability_matrix = zeros(length(participants), length(participants));


        for i = 1:length(participants) 
            for j = 1:length(participants) 
                % Calculate Pearson correlation
                Pearson_Identifiability_matrix(i, j) = corr2(timepoint1{i}, timepoint2{j});
            end
        end

        for i = 1:length(participants) 
            for j = 1:length(participants)
                % Extract the functional connectivity matrices for participants i and j
                A = timepoint1{i};
                B = timepoint2{j};
        
                % Compute the covariance-based identifiability between A and B
                covstatis_Identifiability_matrix(i, j) = (trace(A' * B)) / sqrt(trace(A' * A) * trace(B' * B));
            end
        end
        
        resultStrings = strings(1, length(participants)); % Initialize a string array to store the results

        % For Pearson

        Pearson_Dself_values = zeros(1, length(participants)); % Initialize an array to store Dself values
        
        for p = 1:length(participants)
            % Calculate the differentiability of a participant
            
            % Take the correlation of the participant with itself over time
            Pearson_selfcorr = Pearson_Identifiability_matrix(p, p);
            
            % Take the mean of the correlation of the participant with the others
            Pearson_column_elements = Pearson_Identifiability_matrix(:, p);
            Pearson_column_elements_without_diagonal = Pearson_column_elements([1:p-1, p+1:end]);
            Pearson_meanothercorr = mean(Pearson_column_elements_without_diagonal);
            
            % Calculate the empirical standard deviation of inter-individual features correlations
            Pearson_triangle_identifiability_matrix = tril(Pearson_Identifiability_matrix, -1);
            
            % Calculate the standard deviation of inter-individual feature correlations
            Pearson_sd = std(Pearson_triangle_identifiability_matrix(:));
            
            Pearson_Dself = (Pearson_selfcorr - Pearson_meanothercorr) / Pearson_sd;
            
            % Store the result in the arrays
            resultStrings(p) = sprintf('Participant %d: %f', p, Pearson_Dself);
            Pearson_Dself_values(p) = Pearson_Dself;
        end
        
        % Calculate the mean and median of Dself values
        Pearson_mean_Dself = mean(Pearson_Dself_values);
        Pearson_median_Dself = median(Pearson_Dself_values);
        
        % Display mean and median
        disp(['Pearson_Mean Dself: ' num2str(Pearson_mean_Dself)]);
        disp(['Pearson_Median Dself: ' num2str(Pearson_median_Dself)]);
        
        % For covstatis
        
        covstatis_Dself_values = zeros(1, length(participants)); % Initialize an array to store Dself values
        
        for p = 1:length(participants)
            % Calculate the differentiability of a participant
            
            % Take the correlation of the participant with itself over time
            covstatis_selfcorr = covstatis_Identifiability_matrix(p, p);
            
            % Take the mean of the correlation of the participant with the others
            covstatis_column_elements = covstatis_Identifiability_matrix(:, p);
            covstatis_column_elements_without_diagonal = covstatis_column_elements([1:p-1, p+1:end]);
            covstatis_meanothercorr = mean(covstatis_column_elements_without_diagonal);
            
            % Calculate the empirical standard deviation of inter-individual features correlations
            covstatis_triangle_identifiability_matrix = tril(covstatis_Identifiability_matrix, -1);
            
            % Calculate the standard deviation of inter-individual feature correlations
            covstatis_sd = std(covstatis_triangle_identifiability_matrix(:));
            
            covstatis_Dself = (covstatis_selfcorr - covstatis_meanothercorr) / covstatis_sd;
            
            % Store the result in the arrays
            resultStrings(p) = sprintf('Participant %d: %f', p, covstatis_Dself);
            covstatis_Dself_values(p) = covstatis_Dself;
        end
        
        % Calculate the mean and median of Dself values
        covstatis_mean_Dself = mean(covstatis_Dself_values);
        covstatis_median_Dself = median(covstatis_Dself_values);
        
        % Display mean and median
        disp(['covstatis_Mean Dself: ' num2str(covstatis_mean_Dself)]);
        disp(['covstatis_Median Dself: ' num2str(covstatis_median_Dself)]);       
        
        % Return the result strings
        dasilva = resultStrings;

end

