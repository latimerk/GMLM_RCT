dx = 1e-3;
zz = -10:dx:10;

aa = 1.2;
yy = 0.9;


ff = zeros(numel(zz),1);
gg = zeros(numel(zz),1);

for ii = 1:numel(zz)
    [ff(ii),g] = llFunc([zz(ii);aa], yy);
    gg(ii) = g(1);
end


ff2 = zeros(numel(zz),1);
gg2 = zeros(numel(zz),1);

for ii = 1:numel(zz)
    [ff2(ii),g] = llFunc([aa;zz(ii)], yy);
    gg2(ii) = g(2);
end

gg_est  = (ff(3:end) - ff(1:end-2))./(2*dx);
gg2_est = (ff2(3:end) - ff2(1:end-2))./(2*dx);

figure(1);
clf
subplot(1,2,1);
hold on
plot(zz,gg);
plot(zz(2:end-1), gg_est);


subplot(1,2,2);
hold on
plot(zz,gg2);
plot(zz(2:end-1), gg2_est);