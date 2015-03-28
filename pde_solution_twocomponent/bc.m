function [ pl,ql,pr,qr ] = bc( xl,ul,xr,ur,t,Da,Dr,l,L,m0,K,n,kon_a,koff_a,kon_r,koff_r)
    pl=[0;0;ul(3)];
    pr=[0;0;ur(3)];
    ql=[1;1;0];
    qr=[1;1;0];
end