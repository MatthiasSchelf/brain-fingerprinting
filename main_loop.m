%brainstorm --nogui
gainonly = 0;
datadir='/home/daniele.marinazzo@rc.fhsc.net/Documents/fp_subjects/';
if gainonly
    resdir='/home/daniele.marinazzo@rc.fhsc.net/Documents/fp_results_LF/';
else
    resdir='/home/daniele.marinazzo@rc.fhsc.net/Documents/fp_results/';
end

atlasData = load([datadir 'scout_Desikan-Killiany_68.mat']);
atlas = {atlasData.Scouts.Vertices};

allFiles = dir(fullfile(datadir, '*HEAD*'));
par_0 = 1;
fixed_head = [datadir 'Par' num2str(par_0,'%03d') '_HEAD.mat'];
fixed_data_TX = [datadir 'Par' num2str(par_0,'%03d') '_TX.mat'];

fixed_head_name = ['Par' num2str(par_0,'%03d') '_HEAD.mat'];
allFiles = allFiles(~strcmp({allFiles.name}, fixed_head_name));
participants = {allFiles.name};

fq_tot = [0.5,4,8,13,30; 4,8,13,30,48];
cond_lab = {'matched', 'onehead', 'onedata'};
freq_lab = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
sen_sour_lab = {'sensors','sources'};
functions = {@AEC, 'AEC'; @AEC_noorth, 'AEC_noorth'; @ciPLV, 'ciPLV'};
fs = 200;

parfor isub = 1:length(participants)
    process_subject(isub, participants, datadir, resdir, ...
        fixed_head, fixed_data_TX, atlas, ...
        cond_lab, freq_lab, sen_sour_lab, ...
        fq_tot, functions, fs, gainonly);
end

function process_subject(isub, participants, datadir, resdir, ...
    fixed_head, fixed_data_TX, atlas, ...
    cond_lab, freq_lab, sen_sour_lab, ...
    fq_tot, functions, fs, gainonly)

nconds = length(cond_lab);
nfreqs = length(freq_lab);

for icond = 1%nconds:-1:1
    condition = cond_lab{icond};
    switch condition
        case 'matched'
            head = [datadir participants{isub}];
            for itime = 1:2
                data = [datadir strrep(participants{isub}, 'HEAD', ['T' num2str(itime)])];
                dataData = load(data);
                F = dataData.F(1:306, 201:end);
                F = resample(F', 1, 3)';
                Ws = [1, 40];
                [b, a] = butter(3, Ws / (fs / 2), 'bandpass');
                %filtered = filtfilt(b, a, F');
                %FC = corrcoef(filtered);
                %savefile = [resdir 'corr_wideband_sensors_' condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                %parsave(savefile, FC);

                %sources_avg = getsources(head, data, atlas, gainonly);
                %filtered = filtfilt(b, a, sources_avg);
                %FC = corrcoef(filtered);
                %savefile = [resdir 'corr_wideband_sources_' condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                %parsave(savefile, FC);

                for ifreq = 1:nfreqs
                    functions_local = functions;
                    Ws = fq_tot(:,ifreq)';
                    [b, a] = butter(3, Ws / (fs/2), 'bandpass');
                    %filtered_sources = filtfilt(b, a, sources_avg);
                    filtered_sensors = filtfilt(b, a, F');

                    for i = 1%:size(functions_local, 1)
                        func = functions_local{i,1};
                        func_name = functions_local{i,2};
                        for i_sen_sour = 1%:2
                            if i_sen_sour == 1
                                input_data = filtered_sensors;
                            else
                                input_data = filtered_sources;
                            end
                            if i==1 && i_sen_sour==1
                                continue
                            end
                            FC = func(input_data);
                            savefile = [resdir func_name '_' freq_lab{ifreq} '_' ...
                                sen_sour_lab{i_sen_sour} '_' condition ...
                                '_' num2str(isub,'%03d') '_T' num2str(itime)];
                            parsave(savefile, FC);
                        end
                    end
                end
            end

        case 'onehead'
            head = fixed_head;
            for itime = 1:2
                data = [datadir strrep(participants{isub}, 'HEAD', ['T' num2str(itime)])];
                %sources_avg = getsources(head, data, atlas, gainonly);
                %[b, a] = butter(3, [1, 40] / (fs/2), 'bandpass');
                %filtered = filtfilt(b, a, sources_avg);
                %FC = corrcoef(filtered);
                %savefile = [resdir 'corr_wideband_sources_' condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                %parsave(savefile, FC);

                for ifreq = 1:nfreqs
                    functions_local = functions;
                    Ws = fq_tot(:,ifreq)';
                    [b, a] = butter(3, Ws / (fs/2), 'bandpass');
                    %filtered_sources = filtfilt(b, a, sources_avg);

                    for i = 1%:size(functions_local,1)
                        func = functions_local{i,1};
                        func_name = functions_local{i,2};
                        FC = func(filtered_sources);
                        savefile = [resdir func_name '_' freq_lab{ifreq} '_sources_' ...
                            condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                        parsave(savefile, FC);
                    end
                end
            end

        case 'onedata'
            head = [datadir participants{isub}];
            head_zero = fixed_head;
            for itime = 1:2
                data_zero = strrep(fixed_data_TX, 'X', num2str(itime));
                sources_zero = getsources(head_zero, data_zero, {});
                F = forwardmodel(head, sources_zero);
                data.F = F;
                [b, a] = butter(3, [1, 40] / (fs/2), 'bandpass');
                filtered = filtfilt(b, a, F');
                FC = corrcoef(filtered);
                savefile = [resdir 'corr_wideband_sensors_' condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                %parsave(savefile, FC);

                sources_avg = getsources(head, data, atlas, gainonly);
                filtered = filtfilt(b, a, sources_avg);
                FC = corrcoef(filtered);
                savefile = [resdir 'corr_wideband_sources_' condition '_' num2str(isub,'%03d') '_T' num2str(itime)];
                %parsave(savefile, FC);

                for ifreq = 1:nfreqs
                    functions_local = functions;
                    Ws = fq_tot(:,ifreq)';
                    [b, a] = butter(3, Ws / (fs/2), 'bandpass');
                    filtered_sources = filtfilt(b, a, sources_avg);
                    filtered_sensors = filtfilt(b, a, F');
                    for i = 1%:size(functions_local,1)
                        func = functions_local{i,1};
                        func_name = functions_local{i,2};
                        for i_sen_sour = 1:2
                            if i_sen_sour == 1
                                input_data = filtered_sensors;
                            else
                                input_data = filtered_sources;
                            end
                            if i==1 && i_sen_sour==1
                                continue
                            end
                            FC = func(input_data);
                            savefile = [resdir func_name '_' freq_lab{ifreq} '_' ...
                                sen_sour_lab{i_sen_sour} '_' condition ...
                                '_' num2str(isub,'%03d') '_T' num2str(itime)];
                            parsave(savefile, FC);
                        end
                    end
                end
            end
    end
end
end

function parsave(filename, FC)
save(filename, 'FC');
end
