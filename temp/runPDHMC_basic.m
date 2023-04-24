function [samples,M] = runPDHMC_basic(HMC_settings, yy)
if(nargin < 2)
    yy = 0.1;
end
figNum = 2;

TotalSamples = HMC_settings.nWarmup + HMC_settings.nSamples;
theta = [0;0];

samples.thetas = nan(2, TotalSamples);
samples.thetas(:,1) = theta;
samples.log_p_accept = nan(1, TotalSamples);
samples.log_post = nan(1, TotalSamples);
samples.accepted = nan(1, TotalSamples);
samples.errors = false(1, TotalSamples);
samples.e = nan(2, TotalSamples);


HMC_state.stepSize.e       = HMC_settings.stepSize.e_0;
HMC_state.stepSize.e_bar   = HMC_settings.stepSize.e_0;
HMC_state.stepSize.x_bar_t = 0;
HMC_state.stepSize.x_t     = 0;
HMC_state.stepSize.H_sum   = 0;
HMC_state.steps            = min(HMC_settings.stepSize.maxSteps, ceil(HMC_settings.stepSize.stepL / HMC_state.stepSize.e));

HMC_state.fpi_1 = HMC_settings.fpi_1;
HMC_state.fpi_2 = HMC_settings.fpi_2;

nlpostFunction = @(theta) llFunc(theta, yy);
samples.log_post(1) = -nlpostFunction(theta);
samples.e(:,1)      = HMC_state.stepSize.e;

M = [1;0.01];


for sample_idx = 2:TotalSamples
    
    Gfunc = @(ww,pp) getGstep(ww, M, pp);
    [samples.accepted(sample_idx), samples.errors(sample_idx), theta_new, samples.log_p_accept(sample_idx)] = PDHMCstep_diag(theta, Gfunc, nlpostFunction, HMC_state);
    theta = theta_new;

    samples.thetas(:,sample_idx) = theta;
    samples.log_post(sample_idx) = -nlpostFunction(theta);

    % adjust step size: during warmup
    HMC_state = kgmlm.fittingTools.adjustHMCstepSize(sample_idx, HMC_state, HMC_settings.stepSize, samples.log_p_accept(sample_idx));
    samples.e(:,sample_idx) = [HMC_state.stepSize.e; HMC_state.stepSize.e_bar];

    
    %% updates the covariance matrix of the hyperparameters
    if(ismember(sample_idx, HMC_settings.M_est.samples ) ) 
        start_idx = HMC_settings.M_est.first_sample(HMC_settings.M_est.samples == sample_idx);
        ww = start_idx:sample_idx;
        %diagonal only
        M = (1./var(samples.thetas(:,ww),[],2));
    end

    %% plot updates
    if( mod(sample_idx,500) == 0)
        if(sample_idx == TotalSamples)
            ww = (HMC_settings.nWarmup+1):sample_idx;
        else
            ww = max(2,sample_idx-99):sample_idx;
        end
        
        accept_rate = mean(samples.accepted(ww))*100;
        fprintf("HMC step %d / %d (accept per. = %.1f in last %d steps\n", sample_idx, TotalSamples, accept_rate, numel(ww));
        fprintf("\tcurrent step size = %e, HMC steps = %d, num HMC early rejects = %d\n", HMC_state.stepSize.e, HMC_state.steps, sum(samples.errors, "omitnan"));

        %%
        kgmlm.utils.sfigure(figNum);
        clf;
        NR = 2;
        NC = 2;
        subplot(NR, NC, 1);
        plot(samples.thetas(1,1:sample_idx));
        title("mu");
        subplot(NR, NC, 2);
        plot(samples.thetas(2,1:sample_idx));
        title("log std");
        subplot(NR, NC, 3);
        plot(1:sample_idx,samples.e(1:2,1:sample_idx));
        title("step size");
        subplot(NR, NC, 4);
        plot(1:sample_idx,samples.log_post(1:sample_idx));
        title("log post");
        

        drawnow;
    end

end