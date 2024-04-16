function amico = amico(results, participants)
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
    
    % Store contents of timepoint1 and timepoint2 in matrices to check what
    % they look like
    timepoint1_matrix = cell2mat(timepoint1);
    timepoint2_matrix = cell2mat(timepoint2);

    % Here the identifiability matrix is the same
    % Make identifiability matrix
    Pearson_Identifiability_matrix = zeros(length(participants), length(participants)); 
    covstatis_Identifiability_matrix = zeros(length(participants), length(participants)); 

    % Store Identifiability matrix in the workspace
    assignin('base', 'Pearson_Identifiability_matrix', Pearson_Identifiability_matrix);
    assignin('base', 'covstatis_Identifiability_matrix', covstatis_Identifiability_matrix);
        
    for i = 1:length(participants) 
        for j = 1:length(participants)
            % Inside the identifiability matrix calculation loop
            Pearson_Identifiability_matrix(i, j) = corr(cell2mat(timepoint1(i)), cell2mat(timepoint2(j)));
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
    % title('Pearson_Identifiability_Matrix');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;
    % 
    % % Visualize Identifiability matrix
    % figure;
    % imagesc(covstatis_Identifiability_matrix);
    % colorbar;
    % clim([0, 1]); % Set color limits from 0 to 1
    % title('Covstatis_Identifiability_Matrix');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;

    % Store Identifiability matrix in the workspace
    assignin('base', 'Pearson_Identifiability_matrix', Pearson_Identifiability_matrix);
    assignin('base', 'Covstatis_Identifiability_matrix', covstatis_Identifiability_matrix);
    
    % Decide Iself Iothers and Idiff for Pearson 
        
    Pearson_Iself = mean(diag(Pearson_Identifiability_matrix));
    disp(['Pearson Iself: ', num2str(Pearson_Iself)]);
        
    triangle_identifiability_matrix = tril(Pearson_Identifiability_matrix, -1);
    triangle_identifiability_matrix = nonzeros(triangle_identifiability_matrix);
    Pearson_Iothers = mean(triangle_identifiability_matrix(:));
    disp(['Pearson Iothers: ', num2str(Pearson_Iothers)]);
        
    Pearson_Idiff = (Pearson_Iself-Pearson_Iothers)*100;
    disp(['Pearson Idiff: ', num2str(Pearson_Idiff)]);

    % Decide Iself Iothers and Idiff for covstatis 
        
    covstatis_Iself = mean(diag(covstatis_Identifiability_matrix));
    disp(['Covstatis Iself: ', num2str(covstatis_Iself)]);
        
    covstatis_triangle_identifiability_matrix = tril(covstatis_Identifiability_matrix, -1);
    covstatis_triangle_identifiability_matrix = nonzeros(covstatis_triangle_identifiability_matrix);
    covstatis_Iothers = mean(covstatis_triangle_identifiability_matrix(:));
    disp(['Covstatis Iothers: ', num2str(covstatis_Iothers)]);
        
    covstatis_Idiff = (covstatis_Iself-covstatis_Iothers)*100;
    disp(['Covstatis Idiff: ', num2str(covstatis_Idiff)]);
    
    %Indicate what the result needs to be
    results = ["Pearson_Iself =" Pearson_Iself; "Pearson_Iothers =" Pearson_Iothers; "Pearson_Idiff =" Pearson_Idiff;"covstatis_Iself =" covstatis_Iself; "covstatis_Iothers =" covstatis_Iothers; "covstatis_Idiff =" covstatis_Idiff];
    
    % Store timepoint1_matrix and timepoint2_matrix in the workspace
    assignin('base', 'timepoint1_matrix', timepoint1_matrix);
    assignin('base', 'timepoint2_matrix', timepoint2_matrix);
    
    amico = results;

end



