function [res] = SOR_factor_calculator(lx,ly)
    res = 2/(1+sin(acos((cos(pi/lx)+cos(pi/ly))/2)));
end

