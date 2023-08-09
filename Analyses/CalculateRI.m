function riVal = CalculateRI(responseMatrix)
if sum(responseMatrix(1,:)) == 1 || sum(responseMatrix(2,:)) == 1
    riVal = nan;
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
    riVal = (h+fa-1)/(1-((h-fa)^2));
end

end