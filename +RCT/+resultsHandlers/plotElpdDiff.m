function [] = plotElpdDiff(CC2, CC1, fl, sg, NSEs)
addpath ../dmc_orientation/code/figureMaking/
RCT.utils.setupGMLMpaths

if(nargin < 5 || isempty(NSEs))
    NSEs = 2;
end
if(nargin < 4 || isempty(sg))
    sg = 1;
end
if(ismember(sg,[-1 2]))
    sg = -1;
elseif(ismember(sg,[1]))
    sg = 1;
else
    error("Invalid sign option");
end

if(nargin < 2 || isempty(CC1))
    CC1.X2 = CC2(1).X2*0;
    CC1.Z  = CC2(1).Z*0;
    CC1.tts  = CC2(1).tts;
    CC1.Ts_name  = CC2(1).Ts_name;
end


figure(1);
clf
NR = 2;
NC = 4;

colors = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
NM = numel(CC2);
dd_0 = zeros([size(CC1.X2) NM]);

for mm = 1:NM
    dd_0(:,:,:,mm) =  sg*(CC2(mm).X2 - CC1.X2);
end
dd = dd_0;
if(nargin > 2 && ~isempty(fl))
    for mm = 1:NM
        ff = ones(fl,1);
        for ii = 1:size(dd,3)
            dd(:,:,ii,mm) = conv2(dd_0(:,:,ii,mm),ff,"same");
        end
    end
else
    fl = 1;
end

gg = squeeze(nanmean(dd,2));
gg_se = squeeze(nanstd(dd,[],2)./sqrt(sum(~isnan(dd),2))) * NSEs;

roundInt = 20;
lw_0 = 0.1;
lw = 1;

pps = [1 2 3 5 6];



rrs = false(numel(CC1.tts), 5);
yr = [-0.05 0.05];
for ii = 1:5

    switch CC1.Ts_name(ii)
        case "noise_on"
            rr = CC1.tts >= -50 & CC1.tts <= 250;
        case "stim_on"
            rr = CC1.tts >= -50 & CC1.tts <= 250;
        case "targets_on"
            rr = CC1.tts >= -50 & CC1.tts <= 250;
        case "fix_off"
            rr = CC1.tts >= -50 & CC1.tts <= 250;
        case "saccade_end"
            rr = CC1.tts >= -250 & CC1.tts <= 50;
    end
    m1 = min(0,min(gg(rr,ii,:) - gg_se(rr,ii,:),[],"all","omitnan"));
    m2 = max(0,max(gg(rr,ii,:) + gg_se(rr,ii,:),[],"all","omitnan"));

    m1 = floor(m1*roundInt)./roundInt;
    m2 = ceil(m2*roundInt)./roundInt;
    yr(1) = min(m1, yr(1));
    yr(2) = max(m2, yr(2));
    rrs(:,ii) = rr;
end

for ii = 1:5
    subplot(NR, NC, pps(ii))
    hold on
    rr = rrs(:,ii);
    tts = CC1.tts(rr);

    if(NM == 0)
        plot(tts, dd(rr, :, ii), 'color', ones(1,3)*0.9, 'linewidth', lw_0);
    end

    plot([tts(1) tts(end)], [0 0], 'k:');



    plot([0 0], yr, 'k:');


    for mm = 1:NM
        yy = gg(rr,ii,mm);
        se = gg_se(rr,ii,mm);

        vv = ~isnan(yy);
        DMC.plottingTools.plotLineWithErrorRegion(gca, tts(vv), yy(vv), yy(vv) - se(vv), yy(vv) + se(vv), colors(mm,:), lw);
    end
    ylim(yr);
    xlim([tts(1) tts(end)]);

    xlabel(sprintf("time from %s (ms)", CC1.Ts_name(ii)),"interpreter","none");
    ylabel(sprintf("mean log like diff (per trial per %d ms)", fl));

    hold off
end


A = zeros(numel(CC1.Z), NM);
for mm = 1:NM
    A(:,mm) = sg*(CC2(mm).Z - CC1.Z);
end
subplot(NR, NC, 4);
hist(A)
xlabel("LL diff (A-B)");
ylabel("N trials");
if(NM == 1)
    title(sprintf("mean = %.2f, median = %.2f, std = %.2f, p = %.3e", mean(A), median(A), std(A), signrank(A)));
else
    for mm = 1:NM
        fprintf("mm = %d\tmean = %.2f, median = %.2f, std = %.2f, p = %.3e\n", mm, mean(A(:,mm)), median(A(:,mm)), std(A(:,mm)), signrank(A(:,mm)));
    end
end

subplot(NR, NC, NC+4);
hold on
plot(A)
ff2 = ones(5,1)./5;
cc = conv2(A,ff2,"same");
plot(3:(size(A,1)-2), cc(3:end-2,:), 'k:')
xlabel("Trial");
ylabel("LL diff (A-B)");
hold off


setDefaultAxisProperties(gcf);
