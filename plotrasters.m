close all;

region = "SC";
[data_p, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", region);

P = size(data_p.trials(1).Y.(region),2);

NR = 2;
NC = 3;
NF = ceil(P / (NR * NC));
vv = [data_p.trials(:).targets_vertical];
NT = sum(vv);
Xs = zeros(P, NT);

for pp = 1:P
    ff = ceil(pp / (NR*NC));
    fh = figure(ff);
    sp = mod((pp-1),NR*NC)+1;
    subplot(NR,NC,sp)
    X = cell2mat(arrayfun(@(tt) tt.Y.(region)(:,pp), data_p.trials(vv==1), 'UniformOutput', false)');
    Xs(pp,:) = sum(X);
    [r,c] = find(X>0);
    plot(r,c,'o','markersize', 2, 'markerfacecolor', [0 0 0], 'MarkerEdgeColor', 'none')
    title(sprintf("%d: %s", pp, data_p.NeuronInfo.neuronID(pp)), 'Interpreter', 'none')
    fh.WindowState = 'maximized';
end


NR = 4;
NC = 6;
NF = ceil(P / (NR * NC));


for pp = 1:P
    ff = ceil(pp / (NR*NC));
    fh = figure(ff+100);
    sp = mod((pp-1),NR*NC)+1;
    subplot(NR,NC,sp)
    plot(Xs(pp,:));
    title(sprintf("%d: %s", pp, data_p.NeuronInfo.neuronID(pp)), 'Interpreter', 'none')
    fh.WindowState = 'maximized';
end

figure(200);
clf;
hold on
plot(1:NT, Xs)
plot(1:NT, mean(Xs,1), 'k', 'linewidth', 2);