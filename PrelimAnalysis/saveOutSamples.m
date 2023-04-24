baseFolder = "~/gitCode/GMLM_RCT/Results/GMLM/Coupling_all/Mingus/";

fs = ["SAMPLES_QUARTER123_RCT_GMLM_Mingus_20201211_LIP_targV_Rs_S32_R32_C32_iC16_LC_iC_FEF_SC_run1.mat_part.mat"];

sn = -10000:10:-1;

for ff = 1:numel(fs)
    fname_samples = sprintf("%s/%s", baseFolder, fs(ff));
    if(exist(fname_samples, "file"))
        load(fname_samples, "sample_idx");
    
        sample_idxs = sn + sample_idx;
        [samples] = RCT.resultsHandlers.loadSample(fname_samples, sample_idxs);
    
        save(sprintf("PrelimAnalysis/s%d.mat", ff), "-v7.3", "samples");
    end
end