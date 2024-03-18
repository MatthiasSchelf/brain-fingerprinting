function results = applyPCAsensors(results, num_components)
    all_heads = fieldnames(results);

    % Initialize an empty matrix to store data from all participants
    all_data_matrix = [];

    % Loop through participants
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
                if contains(lower(current_form), 'sensors')
                    % Extract matrices from the nested structure
                    PCA_vector_matrix = end_result;

                    % Concatenate data from all participants horizontally
                    all_data_matrix = [all_data_matrix, PCA_vector_matrix];
                    assignin('base', 'all_data_matrix', all_data_matrix);
                end
            end
        end
    end

    % Center the combined matrix
    mean_all_data_matrix = mean(all_data_matrix);
    centered_all_data_matrix = all_data_matrix - mean_all_data_matrix;

    % Perform PCA on the combined matrix
    [coeff, score, ~] = pca(centered_all_data_matrix);

    % Retain only the specified number of principal components
    coeff = coeff(:, 1:num_components);
    score = score(:, 1:num_components);

    % Reconstruct the data using retained principal components
    reconstructed_data = score * coeff';

    % Reshape the reconstructed data to match the original structure
    start_col = 1;
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
                if contains(lower(current_form), 'sensors')
                    % Determine the number of columns for the current form
                    [~, num_cols] = size(end_result);

                    % Extract and assign the reconstructed data for the current form
                    reconstructed_form_data = reconstructed_data(:, start_col : (start_col + num_cols - 1));
                    results.(current_head).(current_timepoint).([current_form '_PCA']) = reconstructed_form_data;

                    % Update the starting column for the next form
                    start_col = start_col + num_cols;
                end
            end
        end
    end
end
