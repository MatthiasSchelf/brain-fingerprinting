
%% With PCA


function amico = amico_PCA(results_PCA, participants)
    
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
    Identifiability_matrix_PCA = zeros(numParticipants, numParticipants);
    
    % Calculate the correlation for identifiability matrix
    for i = 1:numParticipants
        for j = 1:numParticipants
            Identifiability_matrix_PCA(i, j) = corr(timepoint1_matrix_PCA(i, :)', timepoint2_matrix_PCA(j, :)');
        end
    end

    % Visualize Identifiability matrix
    figure;
    imagesc(Identifiability_matrix_PCA);
    colorbar;
    clim([0, 1]); % Set color limits from 0 to 1
    title('Identifiability Matrix_PCA');
    xlabel('Participant Index');
    ylabel('Participant Index');
    axis square;

    % Store Identifiability matrix in the workspace
    assignin('base', 'Identifiability_matrix_PCA', Identifiability_matrix_PCA);
        
    % Decide Iself Iothers and Idiff
        
    Iself = mean(diag(Identifiability_matrix_PCA));
    disp(['The Iself of this test example is ', num2str(Iself)]);
        
    triangle_identifiability_matrix_PCA = tril(Identifiability_matrix_PCA, -1);
    triangle_identifiability_matrix_PCA = nonzeros(triangle_identifiability_matrix_PCA);
    Iothers = mean(triangle_identifiability_matrix_PCA(:));
    disp(['The Iothers of this test example is ', num2str(Iothers)]);
        
    Idiff = (Iself-Iothers)*100;
    disp(['The Idiff of this test example is ', num2str(Idiff)]);
    
    % Indicate what the result needs to be
    results_PCA = ["Iself =" Iself; "Iothers =" Iothers; "Idiff =" Idiff];

    % Store timepoint1_matrix and timepoint2_matrix in the workspace
    assignin('base', 'timepoint1_matrix_PCA', timepoint1_matrix_PCA);
    assignin('base', 'timepoint2_matrix_PCA', timepoint2_matrix_PCA);
    
    amico = results_PCA;

end