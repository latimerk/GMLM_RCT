% saves data from all areas into one file
% spike counts binned into 1 ms bins
function [data] = processRawData(subjects, sessions)
folders = RCT.dataHandlers.getFolders();

if(nargin < 1 || isempty(subjects))
    subjects = folders.subjects;
end

for ss = 1:numel(subjects)
    subject = subjects(ss);
    [~, subject] = RCT.dataHandlers.checkSubject(folders, subject);
    if(nargin < 2 || isempty(sessions) || (iscell(sessions) && isempty(sessions{ss})))
        sessions_c = 1:numel(folders.sessions.(subject));
    else
        sessions_c = sessions;
    end
    
    for ii = 1:numel(sessions_c)
        [data_raw, session_name] = RCT.dataHandlers.loadRawData(folders, subject, sessions_c(ii));
    
        NL = numel(data_raw); % number of recording locations in file
        % makes some empty values nan for vectorization
        fields_to_fix = ["Align_to_fix_on", "trial_error", "category", "Align_to_cat_stim_on", "Align_to_noise_on", "Align_to_targets_on", "Align_to_fix_off", "Align_to_saccade_end", "tar_location", "targets_vert", "saccade_location"];
        for jj = 1:NL
            for tt = 1:numel(data_raw(jj).Bhv.Trial_info)
                for ff = 1:numel(fields_to_fix)
                    if(isempty(data_raw(jj).Bhv.Trial_info(tt).(fields_to_fix(ff))))
                        data_raw(jj).Bhv.Trial_info(tt).(fields_to_fix(ff)) = nan;
                    end
                end
            end
        end
    
        %% gets trials with consistent timing across all Bhv structs & error codes 0 or 6
        Bhv = data_raw(1).Bhv;
        timings_0 = [Bhv.Trial_info(:).Align_to_fix_on];
        for jj = 2:NL
            timings_c = [data_raw(jj).Bhv.Trial_info(:).Align_to_fix_on];
            timings_0 = timings_0(ismember(timings_0, timings_c));
        end
    
        timings_0 = timings_0(~isnan(timings_0));
        valid_tr_codes = [0 6];
        valid_cats = [-1 1];
    
        timings_c  = [Bhv.Trial_info(:).Align_to_fix_on];
        tr_codes_c = [Bhv.Trial_info(:).trial_error];
        cats_c = [Bhv.Trial_info(:).category];
        vv = ismember(timings_c, timings_0) & ismember(tr_codes_c, valid_tr_codes) & ismember(cats_c, valid_cats);
    
        fs = fieldnames(Bhv);
        for jj = 1:numel(fs)
            if(numel(Bhv.(fs{jj})) == numel(vv))
                Bhv.(fs{jj}) = Bhv.(fs{jj})(vv);
            end
        end
    
        %% log cell info
        data = struct();
        data.bin_size_ms = 1;

        data.session = session_name;
        data.subject = subject; 
        data.NeuronInfo = struct("location", cell(NL,1), "sortingClassification", [], "neuronID", []); 
        NC = nan(NL,1);
        for ll = 1:NL
            data.NeuronInfo(ll).location = data_raw(ll).location;
            NC(ll) = numel(data_raw(ll).Neuro.Neuron);
            data.NeuronInfo(ll).sortingClassification = repmat("", NC(ll),1);
            data.NeuronInfo(ll).neuronID              = repmat("", NC(ll),1);
            for mm = 1:NC(ll)
                data.NeuronInfo(ll).sortingClassification(mm) = lower(string(strtrim(data_raw(ll).Neuro.Neuron(mm).NeuronLabel)));
                data.NeuronInfo(ll).neuronID(mm)              = lower(string(strtrim(data_raw(ll).Neuro.Neuron(mm).NeuronID)));
            end
        end
    
        %% log category info
        dirs = [Bhv.Trial_info(:).direction];
        cats = [Bhv.Trial_info(:).category];
        if(sum(cats == -1) > 0 && sum(cats == 2) > 0)
            error("Unknown category setup.");
        end
        cats(cats == -1) = 2;
        ds = unique(dirs);
        ds = ds(:);
        cs = nan(numel(ds),1);
        for dd = 1:numel(ds)
            uc = unique(cats(dirs == ds(dd)));
            if(numel(uc) ~= 1)
                error("Inconsistent category & direction setup.");
            else
                cs(dd) = uc;
            end
        end
        [cs,or] = sort(cs);
        ds = ds(or);
        ucs = unique(cs);
        for cc = 1:numel(ucs)
            ds_0 = ds(cs == ucs(cc));
            ds_c = mod(ds_0 - 45, 360);
            [~,or] = sort(ds_c);
            ds(cs == ucs(cc)) = ds_0(or);
        end
    
        data.TaskInfo.directions = ds;
        data.TaskInfo.categories = cs;
    
        %% parses each trial
        NT = numel(Bhv.Trial_info);
        data.trials = struct("fix_on", cell(NT,1), "noise_on", [], "targets_on", [], "fixation_time", [], "fix_off", [], "stim_on", [], "saccade_end", [], "saccade_location", [], "target_location", [], "targets_vertical", [], "direction_idx", [], "error_code", [], "Y", []);
        st = 500; % noise on at this time
        TT_max = 3e3 + st;
    
        SAMPLING_RATE = 40;
        BHV_UNITS = 1e3;
    
        for tt = 1:NT
            %% align to noise_on
            t_0 = Bhv.Trial_info(tt).Align_to_noise_on*BHV_UNITS - st;
    
            data.trials(tt).noise_on    = floor(Bhv.Trial_info(tt).Align_to_noise_on   *BHV_UNITS - t_0 + 1);
            data.trials(tt).stim_on     = floor(Bhv.Trial_info(tt).Align_to_cat_stim_on*BHV_UNITS - t_0 + 1);
            data.trials(tt).targets_on  = floor(Bhv.Trial_info(tt).Align_to_targets_on *BHV_UNITS - t_0 + 1);
            data.trials(tt).fix_on      = floor(Bhv.Trial_info(tt).Align_to_fix_on     *BHV_UNITS - t_0 + 1);
            data.trials(tt).fix_off     = floor(Bhv.Trial_info(tt).Align_to_fix_off    *BHV_UNITS - t_0 + 1);
            data.trials(tt).saccade_end = floor(Bhv.Trial_info(tt).Align_to_saccade_end*BHV_UNITS - t_0 + 1);

            data.trials(tt).fixation_time = data.trials(tt).stim_on - 500; %HARDCODED EVENT TIME
    
            data.trials(tt).direction_idx    = find(Bhv.Trial_info(tt).direction == data.TaskInfo.directions,1);
            data.trials(tt).error_code       = Bhv.Trial_info(tt).trial_error;
            data.trials(tt).targets_vertical = Bhv.Trial_info(tt).targets_vert;
            data.trials(tt).target_location  = Bhv.Trial_info(tt).tar_location;
            data.trials(tt).saccade_location = Bhv.Trial_info(tt).saccade_location;
    
            %% get spike times from each recording location
            for ll = 1:NL
                data.trials(tt).Y.(data_raw(ll).location) = zeros(TT_max, NC(ll)); 
    
                for mm = 1:NC(ll)
                    sts = floor(data_raw(ll).Neuro.Neuron(mm).NeuronSpkT/SAMPLING_RATE - t_0 + 1);
                    sts = sts(sts > 0 & sts < TT_max);
                    if(~isempty(sts))
                        data.trials(tt).Y.(data_raw(ll).location)(sts, mm) = 1;
                    end
                end
            end
        end
    
        %% save out session
        fname = RCT.dataHandlers.getProcessedDataFileName(folders, subject, session_name);
        save(fname, "data", "-v7.3");
    
    end
end