function [ ic ] = ic(x,D,l,L,m0,K,n,kon,koff) 
    ic = [kon/koff-0.01*cos(2*pi*x/L); 0];
end

