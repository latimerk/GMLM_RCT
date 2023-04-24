function [HMC_settings] = setupHMCparams()

HMC_settings.stepSize.e_0 = 1e-3;
HMC_settings.stepSize.delta  = 0.8; %0.8
HMC_settings.stepSize.gamma = 0.05;
HMC_settings.stepSize.kappa = 0.75;
HMC_settings.stepSize.t_0   = 10;
HMC_settings.stepSize.mu    = log(10*HMC_settings.stepSize.e_0);
HMC_settings.stepSize.max_step_size = 0.2;


HMC_settings.M_est.first_sample  = [1001];%[201   1501 5001]; %when to estimate cov matrix. At sample=samples(ii), will use first_sample(samples(ii)):sample
HMC_settings.M_est.samples       = [11000];%[1000  4000 20001];

HMC_settings.stepSize.schedule   = [2     24500];


%step size paramters
HMC_settings.stepSize.stepL     = 1.0; %total steps to take is min(maxSteps , ceil(stepL/e))
HMC_settings.stepSize.maxSteps  = 100; %max number of steps per sample


HMC_settings.nWarmup  = 25e3;
HMC_settings.nSamples = 50e3;