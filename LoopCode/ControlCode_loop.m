% run brainstorm first
% Control code loop:
global analyse condition sensors fq_index PCA
%Run brainstorm first


%% Decide the variables:
num_components = 90; % If you want to do PCA decide components
participants = {'Par1', 'Par2', 'Par3', 'Par4', 'Par5', 'Par6', 'Par7', 'Par8', 'Par9', 'Par10', 'Par11', 'Par12', 'Par13', 'Par14', 'Par15', 'Par16', 'Par17', 'Par18', 'Par19', 'Par20', 'Par21', 'Par22', 'Par23', 'Par24', 'Par25', 'Par26', 'Par27', 'Par28', 'Par29', 'Par30', 'Par31', 'Par32', 'Par33', 'Par34', 'Par35', 'Par36', 'Par37', 'Par38', 'Par39', 'Par40', 'Par41', 'Par42', 'Par43', 'Par44', 'Par45', 'Par46', 'Par47', 'Par48', 'Par49', 'Par50', 'Par51', 'Par52', 'Par53', 'Par54', 'Par55', 'Par56', 'Par57', 'Par58', 'Par59', 'Par60', 'Par61', 'Par62', 'Par63', 'Par64', 'Par65', 'Par66', 'Par67', 'Par68', 'Par69', 'Par70', 'Par71', 'Par72', 'Par73', 'Par74', 'Par75', 'Par76', 'Par77', 'Par78', 'Par79', 'Par80', 'Par81', 'Par82', 'Par83', 'Par84', 'Par85', 'Par86', 'Par87', 'Par88', 'Par89', 'Par90', 'Par91', 'Par92', 'Par93', 'Par94', 'Par95', 'Par96', 'Par97', 'Par98', 'Par99', 'Par100', 'Par101', 'Par102', 'Par103', 'Par104', 'Par105', 'Par106', 'Par107', 'Par108', 'Par109', 'Par110', 'Par111', 'Par112', 'Par113', 'Par114', 'Par115', 'Par116', 'Par117', 'Par118', 'Par119', 'Par120'};
fq_tot=[0.5,4,8,13,30;4,8,13,30,48];
fixed_head = 'Par1_Head';
shuffledheads = participants(randperm(length(participants)));
fixed_data_T1 = 'Par1_T1';
fixed_data_T2 = 'Par1_T2';

for PCA = 0:1 %Change if you want to do PCA (0 = no; 1 = yes)
    for sensors = 0:1 % Change if you want to see sensordata: no sensors = 0, sensors = 1
        for condition = 1:4 %Change between Nothing = 0, Matched = 1 or OneHead = 2, Random = 3, OneData = 4
            for analyse = 2:4 %Change to choose the analysis you want 0 = Nothing ,1 = timedomain, 2 = AEC no ort, 3 = ciPLV, 4 = AEC ort
                for fq_index = 1:5 %Decide the frequency for the frequency measures 1=[0.5;4];2=[4;8];3=[8;13];4=[13;30];5=[30;48]
                    fq=fq_tot(:,fq_index);
                    disp([PCA sensors condition analyse fq_index]);
                    %% Run this code once the variables are decided.

                    % Without PCA

                    if PCA == 0
                        if sensors == 0
                            if analyse == 0
                                disp("Put the variables correct.")
                            elseif analyse == 1
                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = conTDMatched(participants);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Matched_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 2
                                    results = conTDOneHead(participants, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneHead_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 3
                                    results = conTDRandom(participants, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Random_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                else
                                    results = conTDOneData(participants, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneData_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            elseif analyse == 2

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECnoortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Matched_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 2
                                    results = AECnoortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneHead_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 3
                                    results = AECnoortRandom(participants, fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Random_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                else
                                    results = AECnoortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneData_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            elseif analyse == 3

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = ciPLVMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Matched_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 2
                                    results = ciPLVOneHead(participants, fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneHead_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 3
                                    results = ciPLVRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Random_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                else
                                    results = ciPLVOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneData_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            else

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Matched_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 2
                                    results = AECortOneHead(participants, fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneHead_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                elseif condition == 3
                                    results = AECortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_Random_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                else
                                    results = AECortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_OneData_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            end

                        else

                            if analyse == 0
                                disp("Put the variables correct.")
                            elseif analyse == 1
                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = conTDMatched(participants);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 2
                                    results = conTDOneHead(participants, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 3
                                    results = conTDRandom(participants, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                else
                                    results = conTDOneData(participants, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            elseif analyse == 2

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECnoortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 2
                                    results = AECnoortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 3
                                    results = AECnoortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                else
                                    results = AECnoortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            elseif analyse == 3

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = ciPLVMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 2
                                    results = ciPLVOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 3
                                    results = ciPLVRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                else
                                    results = ciPLVOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end

                            else

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECortMatched(participants,fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 2
                                    results = AECortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                elseif condition == 3
                                    results = AECortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amicosensors(results,participants);
                                    resultsdasilva = dasilvasensors(results,participants);
                                else
                                    results = AECortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'NonPCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    resultsamico = amico(results,participants);
                                    resultsdasilva = dasilva(results,participants);
                                end
                            end
                        end

                    else

                        %With PCA

                        if sensors == 0
                            if analyse == 0
                                disp("Put the variables correct.")
                            elseif analyse == 1
                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = conTDMatched(participants);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Matched_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = conTDOneHead(participants, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneHead_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = conTDRandom(participants, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Random_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                else
                                    results = conTDOneData(participants, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneData_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            elseif analyse == 2

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECnoortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Matched_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = AECnoortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneHead_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = AECnoortRandom(participants, fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Random_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                else
                                    results = AECnoortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneData_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            elseif analyse == 3

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = ciPLVMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Matched_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = ciPLVOneHead(participants, fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneHead_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = ciPLVRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Random_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                else
                                    results = ciPLVOneData(participants,fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneData_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            else

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Matched_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = AECortOneHead(participants, fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneHead_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = AECortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_Random_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                else
                                    results = AECortOneData(participants,fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_OneData_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            end

                        else

                            if analyse == 0
                                disp("Put the variables correct.")
                            elseif analyse == 1
                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = conTDMatched(participants);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = conTDOneHead(participants, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = conTDRandom(participants, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                else
                                    results = conTDOneData(participants, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_TD.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            elseif analyse == 2

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECnoortMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = AECnoortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = AECnoortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                else
                                    results = AECnoortOneData(participants, fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECnoort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            elseif analyse == 3

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = ciPLVMatched(participants, fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = ciPLVOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condtition == 3
                                    results = ciPLVRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                else
                                    results = ciPLVOneData(participants,fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_ciPLV.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end

                            else

                                if condition == 0
                                    disp("Put the variables correct.")
                                elseif condition == 1
                                    results = AECortMatched(participants,fq);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 2
                                    results = AECortOneHead(participants,fq, fixed_head);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsesnors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                elseif condition == 3
                                    results = AECortRandom(participants,fq, shuffledheads);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCAsensors(results, num_components);
                                    resultsamico = amicosensors_PCA(results_PCA,participants);
                                    resultsdasilva = dasilvasensors_PCA(results_PCA,participants);
                                else
                                    results = AECortOneData(participants,fq, fixed_data_T1, fixed_data_T2);
                                    save(['.\FC\fq_' num2str(fq_index) '_' 'PCA_FC_sensors_AECort.mat'], 'results');
                                    results = trivec(results);
                                    results_PCA = applyPCA(results, num_components);
                                    resultsamico = amico_PCA(results_PCA,participants);
                                    resultsdasilva = dasilva_PCA(results_PCA,participants);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
