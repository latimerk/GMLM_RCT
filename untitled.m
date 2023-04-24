dim_T = 47;
dim_U = 14;
dim_V = 60;

R = 32;

T = randn(dim_T, R);
U = randn(dim_U, R);
V = randn(dim_V, R);

TU = zeros(dim_T*dim_U,R);
for rr = 1:R
    TU(:,rr) = kron(U(:,rr), T(:,rr));
end

TUV = TU*V';