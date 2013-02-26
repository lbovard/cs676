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
        

end