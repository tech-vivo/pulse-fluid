% Will McFadden (wmcfadden)
% returns active fluid velocity solution for parameters q, with different 
% sets of data given in dat
function [ F ] = act_flu_fun(q, dat )
    F = [];
    for n = 1:length(dat)
        F = [F; active_fluid(q,dat{n})];
    end
end