function [dPrm, h, fa] = CalculateDprime(responseMatrix)
%%
if sum(responseMatrix(1,:)) == 1 || sum(responseMatrix(2,:)) == 1
    dPrm = nan;
    h = nan;
    fa = nan;
else
    h = responseMatrix(1,1)/sum(responseMatrix(1,:));
    if h == 1
        h = (sum(responseMatrix(1,:))-1)/sum(responseMatrix(1,:));
    elseif h == 0
        h = 1/sum(responseMatrix(1,:));
    end
    fa = responseMatrix(2,1)/sum(responseMatrix(2,:));
    if fa == 1
        fa = (sum(responseMatrix(2,:))-1)/sum(responseMatrix(2,:));
    elseif fa == 0
        fa = 1/sum(responseMatrix(2,:));
    end
    
    dPrm = norminv(h)-norminv(fa);
end
%
% beta =  exp((norminv(fa)^2 - norminv(h)^2)/2);
%
% c = -.5 * (norminv(h) + (norminv(fa)));
