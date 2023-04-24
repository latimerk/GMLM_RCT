ff = 1;
addpath ~/gitOther/tensor_toolbox/
load(sprintf("PrelimAnalysis/f%d.mat", ff),  "spkHist", "coupling_r", "coupling_1", "coupling_2", "B");

coupling_r = ttm(tensor(coupling_r), B, 1);
coupling_1 = ttm(tensor(coupling_1), B, 1);
coupling_2 = ttm(tensor(coupling_2), B, 1);
spkHist = ttm(tensor(spkHist), B, 1);

S = size(coupling_r,4);
ss = ones(S,1)./S;

coupling_rm = squeeze(ttm(coupling_r, ss', 4));
coupling_1m = squeeze(ttm(coupling_1, ss', 4));
coupling_2m = squeeze(ttm(coupling_2, ss', 4));
spkHist_m = squeeze(ttm(spkHist, ss', 3));