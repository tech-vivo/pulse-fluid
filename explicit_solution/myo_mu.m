% Will McFadden (wmcfadden)
% this function converts concentration to active stress
function mu = myo_mu(q, c )
    mu = q(1)*c.^q(3)./(q(2)^q(3) + c.^q(3));
end