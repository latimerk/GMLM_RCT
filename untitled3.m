bases = RCT.modelBuilder.setupBasis();

NG = numel(metrics.Groups);
NN = 1:3;
figure(1);
clf
subplot(1,NG,1)
[~,aa] = sort(metrics.Groups(1).N(:, sample_idx), "descend");
plot(bases.stimulus.tts_0, bases.stimulus.B * paramStruct2.Groups(1).T{1}(:,aa(NN)))

subplot(1,NG,2)
[~,aa] = sort(metrics.Groups(2).N(:, sample_idx), "descend");
plot(bases.response.tts_0, bases.response.B * paramStruct2.Groups(2).T{1}(:,aa(NN)))

for gg = 3:NG
    subplot(1,NG,3)
    [~,aa] = sort(metrics.Groups(gg).N(:, sample_idx), "descend");
    plot(bases.spkHist.tts_0, bases.spkHist.B * paramStruct2.Groups(gg).T{1}(:,aa(NN)))
end