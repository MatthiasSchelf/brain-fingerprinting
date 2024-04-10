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
        Pearson_Identifiability_matrix_PCA = zeros(numParticipants, numParticipants);
        covstatis_Identifiability_matrix_PCA = zeros(length(participants), length(participants));

        % Calculate the correlation for identifiability matrix
        for i = 1:numParticipants
            for j = 1:numParticipants
                Pearson_Identifiability_matrix_PCA(i, j) = corr(timepoint1_matrix_PCA(i, :)', timepoint2_matrix_PCA(j, :)');
            end
        end
        
        for i = 1:length(participants) 
            for j = 1:length(participants)
                % Extract the functional connectivity matrices for participants i and j
                A = timepoint1_PCA{i};
                B = timepoint2_PCA{j};
        
                % Compute the covariance-based identifiability between A and B
                covstatis_Identifiability_matrix_PCA(i, j) = (trace(A' * B)) / sqrt(trace(A' * A) * trace(B' * B));
            end
        end

        resultStrings = strings(1, length(participants)); % Initialize a string array to store the results
        
        %For Pearson
        for p = 1:length(participants)
            % Calculate the differentiability of a participant
            
            % Take the correlation of the participant with itself over time
            Pearson_selfcorr = Pearson_Identifiability_matrix_PCA(p, p);
            
            % Take the mean of the correlation of the participant with the others
            Pearson_column_elements = Pearson_Identifiability_matrix_PCA(:, p);
            Pearson_column_elements_without_diagonal = Pearson_column_elements([1:p-1, p+1:end]);
            Pearson_meanothercorr = mean(Pearson_column_elements_without_diagonal);
            
            % Calculate the empirical standard deviation of inter-individual features correlations
            Pearson_triangle_identifiability_matrix = tril(Pearson_Identifiability_matrix_PCA, -1);
            
            % Calculate the standard deviation of inter-individual feature correlations
            Pearson_sd = std(Pearson_triangle_identifiability_matrix(:));
            
            Pearson_Dself = (Pearson_selfcorr - Pearson_meanothercorr) / Pearson_sd;
            
        end
        
        % Calculate the mean and median of Dself values
        Pearson_mean_Dself = mean(Pearson_Dself);
        Pearson_median_Dself = median(Pearson_Dself);
        
        % Display mean and median
        disp(['Pearson Mean Dself: ' num2str(Pearson_mean_Dself)]);
        disp(['Pearson Median Dself: ' num2str(Pearson_median_Dself)]);

        %For covstatis
        for p = 1:length(participants)
            % Calculate the differentiability of a participant
            
            % Take the correlation of the participant with itself over time
            covstatis_selfcorr = covstatis_Identifiability_matrix_PCA(p, p);
            
            % Take the mean of the correlation of the participant with the others
            covstatis_column_elements = covstatis_Identifiability_matrix_PCA(:, p);
            covstatis_column_elements_without_diagonal = covstatis_column_elements([1:p-1, p+1:end]);
            covstatis_meanothercorr = mean(covstatis_column_elements_without_diagonal);
            
            % Calculate the empirical standard deviation of inter-individual features correlations
            covstatis_triangle_identifiability_matrix = tril(covstatis_Identifiability_matrix_PCA, -1);
            
            % Calculate the standard deviation of inter-individual feature correlations
            covstatis_sd = std(covstatis_triangle_identifiability_matrix(:));
            
            covstatis_Dself = (covstatis_selfcorr - covstatis_meanothercorr) / covstatis_sd;
            
        end
        
        % Calculate the mean and median of Dself values
        covstatis_mean_Dself = mean(covstatis_Dself);
        covstatis_median_Dself = median(covstatis_Dself);
        
        % Display mean and median
        disp(['covstatis Mean Dself: ' num2str(covstatis_mean_Dself)]);
        disp(['covstatis Median Dself: ' num2str(covstatis_median_Dself)]);
        
        % Return the result strings
        dasilva = resultStrings;
end

