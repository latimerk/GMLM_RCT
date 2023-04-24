%function [CS] = allignTotalLL2(Ts, Ts_name, LL_s, LL, dim_N_ranges)

T_max = 500;
M = numel(dim_N_ranges) - 1;

tts = -T_max:T_max;

K = size(Ts,2);

X2 = nan(T_max*2 + 1, M, K);
X2f = nan(T_max*2 + 1, M, K);
FR = nan(T_max*2 + 1, M, K);
for mm = 1:M

    if(mm == 1 || mod(mm,10) == 0 || mm == M)
        fprintf("Trial %d / %d\n", mm, M);
    end
    tr_length = dim_N_ranges(mm+1) - dim_N_ranges(mm);
    xx = dim_N_ranges(mm):(dim_N_ranges(mm+1)-1);

    for kk = 1:K
        if(~isnan(Ts(mm,kk)))
            tts_c = (1:tr_length) + (T_max + 1 - Ts(mm,kk));
            vv = tts_c >= 1 & tts_c <= numel(tts);

            X2(tts_c(vv), mm, kk) = LL_s(xx(vv));
            X2f(tts_c(vv), mm, kk) = LL2_s(xx(vv));
            FR(tts_c(vv), mm, kk) = FR_s(xx(vv));

        end
    end
end


CS.Z = Z;
CS.tts = tts;
CS.Ts_name = Ts_name;
CS.X2 = X2;
CS.X2f = X2f;
CS.FR = FR;
CS.filterLength = filterLength;