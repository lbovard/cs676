function alphat=interpDelta(alpha,S,St)
    %check if simulated price is outside binomial grid
    %simulated price is lower than smallest binomial grid
    if(St<S(1))
        alphat=alpha(1);
        return;
    end
    if(St>S(end))
        alphat=alpha(end);
        return;
    end
    %otherwise interpolate 
    %matlab profile shows that interp1 is a really slow function
    %easier to just implement linear interpolation oneself
    
    %find the location where S_{i-1} < St < S_{i}
    index=find(logical(S>St),1,'first');
    %Linear slope
    m=(alpha(index)-alpha(index-1))/(S(index)-S(index-1));
    alphat=m*(St-S(index-1))+alpha(index-1);
    
    %alphat=interp1(S,alpha,St);
end
