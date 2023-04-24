% ff = randn(14,1000);
gg = randn(1,1000);
Y = zeros(1, 57); %size(gmlm.X_groups(3).X_local{1},1)

gmlm.setupComputeStructuresHost
% for ii = 1:57
%     Y(ii) = mean(std((gmlm.X_groups(3).X_local{1}(:, (1:14) + (ii-1)*14)*ff).*gg,[],1));
% end

ff2 = randn(size(params.Groups(1).T{1},1),1000);
XX =(GMLMstructure.Groups(1).X_shared{1}'*gmlm.X_groups(1).iX_shared(1).iX)';
Y2 = mean(std(XX*ff2.*gg,[],1))
Y2b = mean(std(XX*ff2,[],1))
Y2c = mean(mean(XX*ff2,1))
Y2d = std(mean(XX*ff2,1))



% figure(1); clf; plot(Y); title(std(Y))