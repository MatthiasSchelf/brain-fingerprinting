function amico = amicosensors(results, participants)
    global analyse condition sensors fq_index PCA
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
    
    % Make identifiability matrix for Pearson
    Pearson_Identifiability_matrix = zeros(length(participants), length(participants));
    covstatis_Identifiability_matrix = zeros(length(participants), length(participants)); 
        
    for i = 1:length(participants) 
        for j = 1:length(participants)
           % Calculate Pearson correlation
           Pearson_Identifiability_matrix(i, j) = corr2(timepoint1{i}, timepoint2{j});
        end
    end
    
    %Make Identifiability matrix for covstatis
    for i = 1:length(participants) 
        for j = 1:length(participants)
            % Extract the functional connectivity matrices for participants i and j
            A = timepoint1{i};
            B = timepoint2{j};
    
            % Compute the covariance-based identifiability between A and B
            covstatis_Identifiability_matrix(i, j) = (trace(A' * B)) / sqrt(trace(A' * A) * trace(B' * B));
        end
    end

    % Save the Identifiability matrix to a MATLAB file
    savename_pearson=['.\IM\Pearson_Identifiability_matrix' '_' num2str(analyse) '_' num2str(condition) '_' num2str(sensors) '_' num2str(fq_index) '_' num2str(PCA)];
    savename_covstatis=['.\IM\covstatis_Identifiability_matrix' '_' num2str(analyse) '_' num2str(condition) '_' num2str(sensors) '_' num2str(fq_index) '_' num2str(PCA)];

    save(savename_pearson, 'Pearson_Identifiability_matrix');
    save(savename_covstatis, 'covstatis_Identifiability_matrix');
    
    % % Visualize Pearson Identifiability matrix
    % figure;
    % imagesc(Pearson_Identifiability_matrix);
    % colorbar;
    % clim([0, 1]); % Set color limits from 0 to 1
    % title('Pearson Identifiability Matrix');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;
    % 
    % % Visualize covstatis Identifiability matrix
    % figure;
    % imagesc(covstatis_Identifiability_matrix);
    % colorbar;
    % clim([0, 1]); % Set color limits from 0 to 1
    % title('covstatis Identifiability Matrix');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;
    
    % Decide Iself Iothers and Idiff for Pearson
        
    Pearson_Iself = mean(diag(Pearson_Identifiability_matrix));
    disp(['Pearson Iself: ', num2str(Pearson_Iself)]);
        
    triangle_identifiability_matrix = tril(Pearson_Identifiability_matrix, -1);
    triangle_identifiability_matrix = nonzeros(triangle_identifiability_matrix);
    Pearson_Iothers = mean(triangle_identifiability_matrix(:));
    disp(['Pearson Iothers: ', num2str(Pearson_Iothers)]);
        
    Pearson_Idiff = (Pearson_Iself-Pearson_Iothers)*100;
    disp(['Pearson Idiff:', num2str(Pearson_Idiff)]);

    % Decide Iself Iothers and Idiff for covstatis
        
    covstatis_Iself = mean(diag(covstatis_Identifiability_matrix));
    disp(['covstatis Iself:', num2str(covstatis_Iself)]);
        
    covstatis_triangle_identifiability_matrix = tril(covstatis_Identifiability_matrix, -1);
    covstatis_triangle_identifiability_matrix = nonzeros(covstatis_triangle_identifiability_matrix);
    covstatis_Iothers = mean(covstatis_triangle_identifiability_matrix(:));
    disp(['covstatis Iothers:', num2str(covstatis_Iothers)]);
        
    covstatis_Idiff = (covstatis_Iself-covstatis_Iothers)*100;
    disp(['covstatis Idiff:', num2str(covstatis_Idiff)]);
    
    %Indicate wha the result needs to be
    results = ["Pearson Iself =" Pearson_Iself; "Iothers =" Pearson_Iothers; "Idiff =" Pearson_Idiff; "covstatis Iself =" covstatis_Iself; "covstatis Iothers =" covstatis_Iothers; "covstatis Idiff =" covstatis_Idiff];
    
    amico = results;

end