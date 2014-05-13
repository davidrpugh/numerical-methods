function [c,k,y,z] = ssrbc(alpha, beta, delta)



check = 0;
k    = (1/alpha*(1/beta - (1-delta)))^(1/(alpha-1));
c     = k^alpha - delta*k;
y     =  k^alpha;
z     = 1;

err(1) =  1/c - beta/c*(alpha*y(+1)/k + 1 - delta);
err(2) =  y - (c + k - (1-delta)*k);
err(3) =  y - (k^alpha)*(z^(1-alpha));


if max(abs(err))>1e-9
    error('fatal (ssrbc1) failed to compute steady state')
end
