function [ fDk , fSk ] = Classify_Dyn_Stat( fS1 , fD1 , p1_1k, p2_1k )
    s_p1 = size(p1_1k);
    for i = 1: s_p1(2):
        idxs = find(fS1==p1_1k(:,i))
        if isempty()
    
    end


end