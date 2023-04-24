%function [LL_a] = allignTotalLL(Ts, LL, LL_prev, LL_curr)

SM = 5;
ff = ones(SM,1)./SM;

T_max = 500; %max(tr_lengths);
M = numel(LL);

tts = -T_max:T_max;

K = size(Ts,2);

NS = size(LL{1},3);
P = size(LL{1},2);

% Y = nan(T_max*2 + 1, M, NS, K);
%W = nan(T_max*2 + 1, M, P, K);
%V = nan(T_max*2 + 1, M, P, K);
X2 = nan(T_max*2 + 1, M, K);


Z = nan(M,1);
for mm = 1:M
    if(mm == 1 || mod(mm,10) == 0 || mm == M)
        fprintf("Trial %d / %d\n", mm, M);
    end
%     LL_tot = squeeze(sum(LL{mm},2));% + LL_prev{mm} - LL_curr{mm};
    LL_tot = kgmlm.utils.logMeanExp(squeeze(sum(LL{mm},2)),2);
    tr_length = size(LL_tot, 1);
    for kk = 1:K
        if(~isnan(Ts(mm,kk)))
            tts_c = (1:tr_length) + (T_max + 1 - Ts(mm,kk));
            vv = tts_c >= 1 & tts_c <= numel(tts);

%             Y_c = reshape(LL_tot, [tr_length 1 NS]);
%             Y(tts_c(vv), mm, :, kk) = Y_c(vv,:);
            X2(tts_c(vv), mm, kk) = LL_tot(vv,:);

            %W(tts_c, mm, :, kk) = reshape(nanmean(SpikeRate{mm}(:,:,:,1),3),  [tr_lengths(mm) 1 P]);
            %V(tts_c, mm, :, kk) = reshape(nanmean(SpikeRate{mm}(:,:,:,2),3),  [tr_lengths(mm) 1 P]);
        end
    end
    Z(mm) = kgmlm.utils.logMeanExp(sum(LL_tot,1),2);
end

%Y2 = squeeze(mean(Y,2,"omitnan")); % mean over trials
%Y3 = squeeze(kgmlm.utils.logMeanExp(Y2,2));


%X2 = squeeze(kgmlm.utils.logMeanExp(Y,3));% mean over samples

%W2 = squeeze(nanmean(W,2));% mean over samples
%V2 = squeeze(nanmean(V,2));% mean over samples

CS.Z = Z;
CS.tts = tts;
%CS.Y3 = Y3;
CS.X2 = X2;
%CS.W2 = W2;
%CS.V2 = V2;
CS.Ts_name = Ts_name;