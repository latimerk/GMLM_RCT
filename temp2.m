bases = RCT.modelBuilder.setupBasis(1, true, true, false);
B = bases.response.B;
k = size(B,2);

C = eye(size(B,1));
D = toeplitz([2;-1;zeros(size(B,1)-2,1)]);
D(1,1) = 1;
D(end,end) = 1;
D = D./2;

figure(1);
clf
hold on
as = 0.1:0.1:0.9;
F = zeros(size(B,1), 1000, numel(as));
for ii = 1:numel(as)
    a = as(ii);

    D2 = toeplitz([2;-a;zeros(size(B,1)-2,1)]);
    D2(1,1) = 1;
    D2(end,end) = 1;
    D2 = D2./2;
    sig2   = inv(B'*((D2 )*B));
    yy = mvnrnd(zeros(k,1), sig2,1000);
    
    F(:,:,ii) = B*yy';
    N = sqrt(sum(F(:,:,ii).^2,1));
    plot(std(F(:,:,ii),[],2))
end