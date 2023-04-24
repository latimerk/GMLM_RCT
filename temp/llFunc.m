function [f,g, aa, bb] = llFunc(theta, y)

nu = 1;
inv_sigma2 = exp(-theta(2)*2);
err = (y - theta(1)).^2;
l_like = -1/2*err.*inv_sigma2 - 1/2*log(2*pi) - theta(2);
lp_mu = -1/2*theta(1).^2 - 1/2*log(2*pi);

if(nargout > 1)
    [lp_var, dp] = DMC.priors.halfTPrior(theta(2), nu);

    g = [(y-theta(1))*inv_sigma2 - theta(1);
         dp + err.*inv_sigma2 - 1];
    g = -g;
else
    lp_var = DMC.priors.halfTPrior(theta(2), nu);
end
f = lp_var + lp_mu + l_like;

aa = [];
bb = [];

f = -f;