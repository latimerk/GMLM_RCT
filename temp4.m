
ff = randn(14,1000);
gg = randn(1,1000);

gmlm.setupComputeStructuresHost;

Y1 = zeros(57,1);
Y2 = zeros(57,1);
for ii = 1:57
xx = gmlm.X_groups(3).X_local{1}(:,(1:14) + (ii-1)*14);

yy = zscore(xx)./sqrt(14);



Y1(ii)= mean(std((xx*ff).*gg,[],1));
Y2(ii)= mean(std((yy*ff).*gg,[],1));
end