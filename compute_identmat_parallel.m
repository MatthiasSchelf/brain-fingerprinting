resdir='/home/daniele.marinazzo@rc.fhsc.net/Documents/fp_results/';

logfile = fullfile(resdir, 'missing_files_log.txt');
fid = fopen(logfile, 'w');

% Gather file names and extract unique prefixes
filelist = dir(fullfile(resdir, '*.mat'));
file_names = {filelist.name};
 
% Extract unique prefixes from filenames
prefixes = {};
for i = 1:length(file_names)
    % Match any characters before subject ID and timepoint pattern
    tokens = regexp(file_names{i}, '^(.*?)_\d{3}_T[12]\.mat$', 'tokens', 'once');
    if ~isempty(tokens)
        prefixes{end+1} = tokens{1};
    end
end
prefixes = unique(prefixes);

fprintf('Found %d unique prefixes for analysis\n', length(prefixes));

% Function to extract upper triangular part (excluding diagonal)
get_upper_tri = @(mat) mat(triu(true(size(mat)), 1));

% Number of subjects
n_subj = 119;

% Initialize parallel pool if not already open
if isempty(gcp('nocreate'))
    parpool('local');
end

% Initialize structure to collect all missing files
all_missing_files = cell(length(prefixes), 1);

% Process each prefix in parallel
parfor i = 1:length(prefixes)
    prefix = prefixes{i};
    fprintf('Processing prefix: %s (%d/%d)\n', prefix, i, length(prefixes));
    
    % Initialize identifiability matrix for this prefix
    ident_mat = nan(n_subj, n_subj);
    
    % Track missing files for this prefix
    missing_files = {};
    
    % Check for expected subject files (1-119)
    for subj = 1:n_subj
        % Check both T1 and T2 files
        t1_filename = sprintf('%s_%03d_T1.mat', prefix, subj);
        t2_filename = sprintf('%s_%03d_T2.mat', prefix, subj);
        
        t1_exists = exist(fullfile(resdir, t1_filename), 'file') == 2;
        t2_exists = exist(fullfile(resdir, t2_filename), 'file') == 2;
        
        % Log if one exists but not the other
        if t1_exists && ~t2_exists
            missing_files{end+1} = t2_filename;
        elseif ~t1_exists && t2_exists
            missing_files{end+1} = t1_filename;
        end
    end
    
    % Pre-load all subject data for this prefix
    vecs_t1 = cell(n_subj, 1);
    vecs_t2 = cell(n_subj, 1);
    valid_t1 = false(n_subj, 1);
    valid_t2 = false(n_subj, 1);
    
    % Load all T1 and T2 data for all subjects
    for subj = 1:n_subj
        % Load T1 data
        t1_filename = sprintf('%s_%03d_T1.mat', prefix, subj);
        t1_filepath = fullfile(resdir, t1_filename);
        if exist(t1_filepath, 'file') == 2
            try
                t1_data = load(t1_filepath, 'FC');
                if isfield(t1_data, 'FC')
                    vecs_t1{subj} = get_upper_tri(t1_data.FC);
                    valid_t1(subj) = true;
                end
            catch
                % Don't log this as a missing file according to your definition
                % Just mark as invalid
            end
        end
        
        % Load T2 data
        t2_filename = sprintf('%s_%03d_T2.mat', prefix, subj);
        t2_filepath = fullfile(resdir, t2_filename);
        if exist(t2_filepath, 'file') == 2
            try
                t2_data = load(t2_filepath, 'FC');
                if isfield(t2_data, 'FC')
                    vecs_t2{subj} = get_upper_tri(t2_data.FC);
                    valid_t2(subj) = true;
                end
            catch
                % Don't log this as a missing file according to your definition
                % Just mark as invalid
            end
        end
    end
    
    % Compute identifiability matrix - correlate all combinations of subjects
    for subj1 = 1:n_subj
        if ~valid_t1(subj1)
            continue;  % Skip if subject doesn't have valid T1 data
        end
        
        vec1 = vecs_t1{subj1};
        
        for subj2 = 1:n_subj
            if ~valid_t2(subj2)
                continue;  % Skip if subject doesn't have valid T2 data
            end
            
            vec2 = vecs_t2{subj2};
            
            % Calculate Pearson correlation between FC vectors
            r = corr(vec1, vec2, 'type', 'Pearson');
            ident_mat(subj1, subj2) = r;
        end
    end
    
    % Save identifiability matrix for this prefix
    outname = sprintf('identmat_%s.mat', prefix);
    parsave(fullfile(resdir, outname), ident_mat);
    
    % Save missing files information
    all_missing_files{i} = missing_files;
end

% Write all missing files to the log outside of parfor
for i = 1:length(all_missing_files)
    missing_files = all_missing_files{i};
    if ~isempty(missing_files)
        fprintf(fid, '=== Missing files for prefix: %s ===\n', prefixes{i});
        for j = 1:length(missing_files)
            fprintf(fid, '%s\n', missing_files{j});
        end
        fprintf(fid, '\n');
    end
end

fclose(fid);
fprintf('Processing complete. Check the log file for any missing files.\n');

% Helper function to save from inside parfor
function parsave(filepath, ident_mat)
    save(filepath, 'ident_mat');
end