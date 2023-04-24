
% NS        = [10   10   20];
% stepSizes = [1e-8 1e-5 1e-3];
NS        = [50];
stepSizes = [1e-3];
nlpostFunction = @(ww) gmlm.vectorizedNLPost_func(ww, params, opts, res);

for ss = 1:numel(stepSizes)
    fprintf("step size %d / %d (%e)\n", ss, numel(stepSizes), stepSizes(ss));
    HMC_state.stepSize.e = stepSizes(ss);
    HMC_state.steps = 100;
    M = gmlm.vectorizeParams(RCT.modelBuilder.getMInit(params, modelSetup, gmlm.dim_M), opts);
    w_init = gmlm.vectorizeParams(params, opts);
    w_c = w_init;
    
    accepted = false(NS(ss),1);
    for ii = 1:NS(ss)
        fprintf("sample %d / %d", ii, NS(ss));
        [accepted(ii), ~, w_c] = kgmlm.fittingTools.HMCstep_diag(w_c,  M, nlpostFunction, HMC_state);
        fprintf("  (accepted = %d)\n", accepted(ii));
    end
    fprintf("Fraction accepted = %.2f\n", mean(accepted));
    w_init = w_c;
    params = gmlm.devectorizeParams(w_c, params, opts);
end