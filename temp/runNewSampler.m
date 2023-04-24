addpath ../../GMLM/example; addpath ../../GMLM_RCT/; addpath ../../GMLM/;

HMC_settings = setupHMCparams();
HMC_settings.fpi_1 = 10;
HMC_settings.fpi_2 = 10;

[samples,M] = runPDHMC_basic(HMC_settings);

%%
figure(11);
clf
ss = samples.thetas(:,25001:end)';
subplot(2,3,1)
plot(ss(:,1), (ss(:,2)), 'o','markersize', 1);
subplot(2,3,2)
plot(ss(:,1), exp(ss(:,2)), 'o','markersize', 1);


subplot(2,3,4)
plot(ss(:,1));
subplot(2,3,5)
plot(ss(:,2));

subplot(2,3,6)
hold on
T = 200;
plot(-T:T, xcorr(zscore(ss(:,1)), T, 'coeff'));
plot(-T:T, xcorr(zscore(ss(:,2)), T, 'coeff'));