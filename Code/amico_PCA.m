
%% With PCA


function amico = amico_PCA(results_PCA, participants)
    global analyse condition sensors fq_index PCA

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
    
                % Check if the current form is based on sensors and contains 'PCA_'
                if ~contains(lower(current_form), 'sensors') && contains(lower(current_form), '_pca')
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

    % Convert cell arrays to matrices for timepoint1 and timepoint2
    timepoint1_matrix_PCA = cell2mat(timepoint1_PCA);
    timepoint2_matrix_PCA = cell2mat(timepoint2_PCA);
    
    % Create an identifiability matrix based on the number of participants
    numParticipants = length(participants);
    Pearson_Identifiability_matrix = zeros(numParticipants, numParticipants); 
    covstatis_Identifiability_matrix = zeros(length(participants), length(participants)); 
    
    % Calculate the correlation for identifiability matrix
    for i = 1:numParticipants
        for j = 1:numParticipants
            Pearson_Identifiability_matrix(i, j) = corr(timepoint1_matrix_PCA(i, :)', timepoint2_matrix_PCA(j, :)');
        end
    end

    for i = 1:length(participants) 
        for j = 1:length(participants)
            % Extract the functional connectivity matrices for participants i and j
            A = timepoint1_PCA{i};
            B = timepoint2_PCA{j};
    
            % Compute the covariance-based identifiability between A and B
            covstatis_Identifiability_matrix(i, j) = (trace(A' * B)) / sqrt(trace(A' * A) * trace(B' * B));
        end
    end

    % Save the Identifiability matrix to a MATLAB file
    savename_pearson=['.\IM\Pearson_Identifiability_matrix' '_' num2str(analyse) '_' num2str(condition) '_' num2str(sensors) '_' num2str(fq_index) '_' num2str(PCA)];
    savename_covstatis=['.\IM\covstatis_Identifiability_matrix' '_' num2str(analyse) '_' num2str(condition) '_' num2str(sensors) '_' num2str(fq_index) '_' num2str(PCA)];

    save(savename_pearson, 'Pearson_Identifiability_matrix');
    save(savename_covstatis, 'covstatis_Identifiability_matrix');

    % % Visualize Identifiability matrix
    % figure;
    % imagesc(Pearson_Identifiability_matrix);
    % colorbar;
    % clim([0, 1]); % Set color limits from 0 to 1
    % title('Pearson_Identifiability Matrix_PCA');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;
    % 
    % % Visualize Identifiability matrix
    % figure;
    % imagesc(covstatis_Identifiability_matrix);
    % colorbar;
    % clim([0, 1]); % Set color limits from 0 to 1
    % title('covstatis_Identifiability Matrix_PCA');
    % xlabel('Participant Index');
    % ylabel('Participant Index');
    % axis square;

    % Store Identifiability matrix in the workspace
    assignin('base', 'Pearson_Identifiability_matrix_PCA', Pearson_Identifiability_matrix);
    assignin('base', 'covstatis_Identifiability_matrix_PCA', covstatis_Identifiability_matrix);
        
    % Decide Iself Iothers and Idiff for Pearsons
        
    Pearson_Iself = mean(diag(Pearson_Identifiability_matrix));
    disp(['Pearson Iself: ', num2str(Pearson_Iself)]);
        
    Pearson_triangle_identifiability_matrix_PCA = tril(Pearson_Identifiability_matrix, -1);
    Pearson_triangle_identifiability_matrix_PCA = nonzeros(Pearson_triangle_identifiability_matrix_PCA);
    Pearson_Iothers = mean(Pearson_triangle_identifiability_matrix_PCA(:));
    disp(['Pearson Iother ', num2str(Pearson_Iothers)]);
        
    Pearson_Idiff = (Pearson_Iself-Pearson_Iothers)*100;
    disp(['Pearson Idiff ', num2str(Pearson_Idiff)]);
    
    % Decide Iself Iothers and Idiff for covstatis
        
    covstatis_Iself = mean(diag(covstatis_Identifiability_matrix));
    disp(['covstatis Iself', num2str(covstatis_Iself)]);
        
    covstatis_triangle_identifiability_matrix_PCA = tril(covstatis_Identifiability_matrix, -1);
    covstatis_triangle_identifiability_matrix_PCA = nonzeros(covstatis_triangle_identifiability_matrix_PCA);
    covstatis_Iothers = mean(covstatis_triangle_identifiability_matrix_PCA(:));
    disp(['covstatis Iothers ', num2str(covstatis_Iothers)]);
        
    covstatis_Idiff = (covstatis_Iself-covstatis_Iothers)*100;
    disp(['covstatis Idiff ', num2str(covstatis_Idiff)]);    
    % Indicate what the result needs to be
    results_PCA = ["Pearson Iself =" Pearson_Iself; "Pearson Iothers =" Pearson_Iothers; "Pearson Idiff =" Pearson_Idiff; "covstatis Iself =" covstatis_Iself; "cvostatis Iothers=" covstatis_Iothers; "covstatis Idiff=" covstatis_Idiff];

    % Store timepoint1_matrix and timepoint2_matrix in the workspace
    assignin('base', 'timepoint1_matrix_PCA', timepoint1_matrix_PCA);
    assignin('base', 'timepoint2_matrix_PCA', timepoint2_matrix_PCA);
    
    amico = results_PCA;

end