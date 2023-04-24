locs = ["LIP" "SC" "FEF"];
if(~exist("trials", "var"))
    [data, session_name] = RCT.dataHandlers.loadData("Mingus", "20201211", locs);
    
    vv = [data.trials(:).targets_vertical];
    trials = data.trials(vv ==1);
    clear data;
end

events = ["noise_on", "noise_on", "targets_on", "fix_off", "saccade_end"];
tts = -500:800;

locs = ["LIP"];
events = [ "saccade_end"];
tts = -200:60;

T = numel(tts);
NT = numel(trials);

ff = ones(5,1);

for jj = 1:numel(locs)
    figure(jj)
    clf

    P = size(trials(1).Y.(locs(jj)),2);
    for bb = 1:numel(events)
        ets = [trials(:).(events(bb))];

        Y = nan(T, P, NT);
        for mm = 1:NT
            if(~isnan(ets(mm)))
                Y(:,:, mm) = trials(mm).Y.(locs(jj))(ets(mm) + tts,:);
            end
        end
        Z = conv2(nanmean(Y,[2 3]),ff,"same");

        mr = [0 ceil(max(Z,[],"all")*50)/50];
        subplot(2,3,bb)
        hold on
        plot([tts(1) tts(end)], [0 0], 'k:');
        plot([0 0], mr, 'k:');
        plot(tts, Z);
        ylim(mr);
        xlabel(events(bb));
        if(bb == 1)
            title(locs(jj))
        end
        hold off
    end
end