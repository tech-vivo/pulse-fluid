function [ ic ] = ic(x,Da,Dr,l,L,m0,K,n,kon_a,koff_a,kon_r,koff_r) 
    ic = [kon_a/koff_a-0.01*cos(2*pi*x/L);kon_r/koff_r-0.01*cos(2*pi*x/L); 0];
end

