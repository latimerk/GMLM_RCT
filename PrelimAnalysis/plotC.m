pp = (1:5) + 25; % destination cells
cc = coupling_2;

NC = 9;
NR = ceil(size(cc,2)/NC);
figure(1);
clf;

is_local = all(double(cc(:,1,1,1)) == 0);
pc = [5 50 95];

colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

for nn = 1:size(cc,2) %over all source cells
    subplot(NR, NC, nn);
    hold on;
    for pp_idx = 1:numel(pp)
        if(is_local && nn == pp(pp_idx))
            yy = prctile(double(squeeze(spkHist(:,pp(pp_idx),:))),pc,2);
            if(numel(pp) == 1)
                title('spk hist');
            end
        else
            yy = prctile(double(squeeze(cc(:,nn,pp(pp_idx),:))),pc,2);
        end
        DMC.plottingTools.plotLineWithErrorRegion(gca, 1:400, yy(:,2), yy(:,1), yy(:,3), colors(mod(pp_idx-1, size(colors,1))+1 ,:));
    end
    plot([0 400], [0 0], 'k--');
end