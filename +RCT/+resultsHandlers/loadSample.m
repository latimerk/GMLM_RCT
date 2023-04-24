function [samples] = loadSample(fname_samples, sample_idx, loadUnScaled)

samples.sample_idx = sample_idx;

if(isstruct(fname_samples) && isfield(fname_samples, "W"))
    samples = fname_samples;
    saveDuring = false;
elseif(isstring(fname_samples) || ischar(fname_samples))
    load(fname_samples, "samples_file_format", "HMC_settings");%"fname_samples_dat"
    if(~isempty(samples_file_format))
        saveDuring = true;
        fname_samples_dat = HMC_settings.samplesFile;
    else
        load(fname_samples, "samples");
        saveDuring = false;
    end
else
    error("fname_samples must be full samples struct or file name to samples .mat file.");
end

if(nargin < 3 || isempty(loadUnScaled))
    loadUnScaled = false;
end

if(saveDuring)
    %% load from binary file
    samples_file = memmapfile(fname_samples_dat,...
                   "Format",samples_file_format, ...
                   "Writable", false);

    if(loadUnScaled && isfield(samples_file.Data, "W_scaled") && ~isempty(samples_file.Data.W_scaled))
        samples.W = samples_file.Data.W_scaled(:,sample_idx);
    else
        samples.W = samples_file.Data.W(:,sample_idx);
    end

    if(loadUnScaled && isfield(samples_file.Data, "B_scaled") && ~isempty(samples_file.Data.B_scaled))
        samples.B = samples_file.Data.B_scaled(:,:,sample_idx);
    elseif(isfield(samples_file.Data, "B") && ~isempty(samples_file.Data.B))
        samples.B = samples_file.Data.B(:,:,sample_idx);
    else
        samples.B = [];
    end
    if(isfield(samples_file.Data, "H") && ~isempty(samples_file.Data.H))
        samples.H = samples_file.Data.H(:,sample_idx);
    else
        samples.H = [];
    end
    if(isfield(samples_file.Data, "H_gibbs") && ~isempty(samples_file.Data.H_gibbs))
        samples.H_gibbs = samples_file.Data.H(:,sample_idx);
    else
        samples.H_gibbs = [];
    end

    jj = 1;
    while(isfield(samples_file.Data, sprintf("G%d_V", jj)))
        if(loadUnScaled && isfield(samples_file.Data, sprintf("G%d_V_scaled", jj)) && ~isempty(samples_file.Data.(sprintf("G%d_V_scaled", jj))))
            samples.Groups(jj).V = samples_file.Data.(sprintf("G%d_V_scaled", jj))(:,:,sample_idx);
        else
            samples.Groups(jj).V = samples_file.Data.(sprintf("G%d_V", jj))(:,:,sample_idx);
        end

        ss = 1;
        while(isfield(samples_file.Data, sprintf("G%d_T_%d", jj, ss)))
            if(loadUnScaled && isfield(samples_file.Data, sprintf("G%d_T_%d_scaled", jj, ss)) && ~isempty(samples_file.Data.(sprintf("G%d_T_%d_scaled", jj, ss))))
                samples.Groups(jj).T{ss} = samples_file.Data.(sprintf("G%d_T_%d_scaled", jj, ss))(:,:,sample_idx);
            else
                samples.Groups(jj).T{ss} = samples_file.Data.(sprintf("G%d_T_%d", jj, ss))(:,:,sample_idx);
            end
            ss = ss + 1;
        end

        if(isfield(samples_file.Data, sprintf("G%d_H", jj)) && ~isempty(samples_file.Data.(sprintf("G%d_H", jj))))
            samples.Groups(jj).H = samples_file.Data.(sprintf("G%d_H", jj))(:,sample_idx);
        end
        if(isfield(samples_file.Data, sprintf("G%d_H_gibbs", jj)) && ~isempty(samples_file.Data.(sprintf("G%d_H_gibbs", jj))))
            samples.Groups(jj).H_gibbs(:) = samples_file.Data.(sprintf("G%d_H_gibbs", jj))(:,sample_idx);
        else
            samples.Groups(jj).H_gibbs = [];
        end

        jj = jj + 1;
    end
else
    %% load from struct
    if(loadUnScaled && isfield(samples, "W_scaled") && ~isempty(samples.W_scaled))
        samples.W = samples.W_scaled(:,sample_idx);
    else
        samples.W = samples.W(:,sample_idx);
    end

    if(loadUnScaled && isfield(samples, "B_scaled") && ~isempty(samples.B_scaled))
        samples.B = samples.B_scaled(:,:,sample_idx);
    else
        samples.B = samples.B(:,:,sample_idx);
    end
    if(isfield(samples, "H") && ~isempty(samples.H))
        samples.H = samples.H(:,sample_idx);
    else
        samples.H = [];
    end
    if(isfield(samples, "H_gibbs") && ~isempty(samples.H_gibbs))
        samples.H_gibbs = samples.H_gibbs(:,sample_idx);
    else
        samples.H_gibbs = [];
    end

    for jj = 1:numel(samples.Groups)
        if(loadUnScaled && isfield(samples.Groups(jj), "V_scaled") && ~isempty(samples.Groups(jj).V_scaled))
            samples.Groups(jj).V = samples.Groups(jj).V_scaled(:,:,sample_idx);
        else
            samples.Groups(jj).V = samples.Groups(jj).V(:,:,sample_idx);
        end

        for ss = 1:numel(samples.Groups.T)
            if(loadUnScaled && isfield(samples.Groups(jj), "T_scaled") && numel(samples.Groups(jj).T_scaled) >= ss &&  ~isempty(samples.Groups(jj).T_scaled{ss}))
                samples.Groups(jj).T{ss} = samples.Groups(jj).T_scaled{ss}(:,:,sample_idx);
            else
                samples.Groups(jj).T{ss} = samples.Groups(jj).T{ss}(:,:,sample_idx);
            end
        end

        if(isfield(samples.Groups(jj), "H") && ~isempty(samples.Groups(jj).H))
            samples.Groups(jj).H = samples.Groups(jj).H(:,sample_idx);
        else
            samples.Groups(jj).H = [];
        end
        if(isfield(samples.Groups(jj), "H_gibbs") && ~isempty(samples.Groups(jj).H_gibbs))
            samples.Groups(jj).H_gibbs = samples.Groups(jj).H_gibbs(:,sample_idx);
        else
            samples.Groups(jj).H_gibbs = [];
        end
    end
end


end