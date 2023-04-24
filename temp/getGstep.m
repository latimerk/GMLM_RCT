function [g, p_ham] = getGstep(theta, G_hat, p)

G_hat = [0.1;0.1];

var_inv = exp(-2*theta(2));

g = [max(G_hat(1), var_inv);
     G_hat(2)];

if(nargout > 1 && G_hat(1) < var_inv)
    var_0 = exp(2*theta(2));
    p_dG_inv_p  = p(1)^2 * var_0 * 2;
    tr_G_inv_dG = -2;

    p_ham = -1/2*tr_G_inv_dG - 1/2*p_dG_inv_p;
    p_ham = [0; p_ham];
else
    p_ham = [0;0];
end