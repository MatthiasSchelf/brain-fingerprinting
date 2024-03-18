function amico = amico(results, participants)
    
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
    Identifiability_matrix = zeros(length(participants), length(participants)); 

    % Store Identifiability matrix in the workspace
    assignin('base', 'Identifiability_matrix', Identifiability_matrix);
        
    for i = 1:length(participants) 
        for j = 1:length(participants)
            % Inside the identifiability matrix calculation loop
            Identifiability_matrix(i, j) = corr(cell2mat(timepoint1(i)), cell2mat(timepoint2(j)));
        end
    end
    
    % Visualize Identifiability matrix
    figure;
    imagesc(Identifiability_matrix);
    colorbar;
    clim([0, 1]); % Set color limits from 0 to 1
    title('Identifiability Matrix');
    xlabel('Participant Index');
    ylabel('Participant Index');
    axis square;


    % Store Identifiability matrix in the workspace
    assignin('base', 'Identifiability_matrix', Identifiability_matrix);
    
    % Decide Iself Iothers and Idiff
        
    Iself = mean(diag(Identifiability_matrix));
    disp(['The Iself of this test example is ', num2str(Iself)]);
        
    triangle_identifiability_matrix = tril(Identifiability_matrix, -1);
    triangle_identifiability_matrix = nonzeros(triangle_identifiability_matrix);
    Iothers = mean(triangle_identifiability_matrix(:));
    disp(['The Iothers of this test example is ', num2str(Iothers)]);
        
    Idiff = (Iself-Iothers)*100;
    disp(['The Idiff of this test example is ', num2str(Idiff)]);
    
    %Indicate what the result needs to be
    results = ["Iself =" Iself; "Iothers =" Iothers; "Idiff =" Idiff];
    
    % Store timepoint1_matrix and timepoint2_matrix in the workspace
    assignin('base', 'timepoint1_matrix', timepoint1_matrix);
    assignin('base', 'timepoint2_matrix', timepoint2_matrix);
    
    amico = results;

end




