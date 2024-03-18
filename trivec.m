function results = trivec(results)
    
    all_heads = fieldnames(results);
    
    for h = 1:length(all_heads)
        current_head = all_heads{h};
        all_timepoints = fieldnames(results.(current_head));
    
        for t = 1:length(all_timepoints)
            current_timepoint = all_timepoints{t};
            all_forms = fieldnames(results.(current_head).(current_timepoint));
    
            for f = 1:length(all_forms)
                current_form = all_forms{f};
                con_matrix = results.(current_head).(current_timepoint).(current_form);
    
                % Take the lower triangle under the diagonal.
                
                triangle = tril(con_matrix, -1);
                
                % Now we vectorize the lower triangle
                
                triangle = triangle(:);
                
                %Indicate what needs to be shown 
                
                trivec = triangle;
    
                % Overwrite the current matrix with the lower triangle vectorized version
                results.(current_head).(current_timepoint).(current_form) = trivec;
            end
        end
    end
    
end