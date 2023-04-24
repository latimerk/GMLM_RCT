function [baselineRates] = getMeanRatesFromBaseline(data, modelSetup_0, sample_idxs)

modelSetup = modelSetup_0;

%% build gmlm for all trials
modelSetup.fitQuarter = 0;
modelSetup.fitHalf = 0;
[GMLMstructure, trials, ~, ~, trialsUsed, ~, ~, timeWindows] = RCT.modelBuilder.constructGMLMdata14(data, modelSetup);

%% get samples
[~, ~, fname_samples] = RCT.modelFitting.getModelFitFileName(data.subject, data.session, modelSetup_0);
samples = RCT.resultsHandlers.loadSample(fname_samples, sample_idxs);

M_max = 100;
M = numel(trials);
M_0 = numel(data.trials);
baselineRates = struct("Y", cell(M_0,1));%, "lambda", cell(M_0,1));

for block = 1:ceil(M/M_max)
    %% build gmlm for subset of trials
    trs_block = (1:M_max) + (block-1)*M_max;
    trs_block = trs_block(trs_block <= M);
    fprintf("Trials %d - %d / %d.\n", trs_block(1), trs_block(end), M);

    gmlm = kgmlm.GMLM(GMLMstructure, trials(trs_block), modelSetup.delta_t_sec, modelSetup.llType);
    
    %% get mean rate
    [FR_s, dim_N_ranges] = RCT.resultsHandlers.evaluateMeanFiringRate(gmlm, samples, modelSetup); %, lambda_s
    
    %% save mean rates
    for mm_i = 1:numel(trs_block)
        mm = trs_block(mm_i);
        ff = trialsUsed(mm);
        baselineRates(ff).Y                    = nan(size(data.trials(ff).Y.(modelSetup.location)));
        baselineRates(ff).Y(timeWindows{mm},:) = FR_s(dim_N_ranges(mm_i):(dim_N_ranges(mm_i+1)-1),:);
%         baselineRates(ff).lambda                    = nan(size(data.trials(ff).Y.(modelSetup.location)));
%         baselineRates(ff).lambda(timeWindows{mm},:) = lambda_s(dim_N_ranges(mm_i):(dim_N_ranges(mm_i+1)-1),:);
    end
end