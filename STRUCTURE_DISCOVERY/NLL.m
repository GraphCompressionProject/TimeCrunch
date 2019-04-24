 function [ l ] = NLL( incl, excl, sub )
        if sub == 0
            l = -log2(excl / (incl + excl));
        elseif sub == 1
            l = -log2(incl / (incl + excl));
        else
            err = 'error... Can only compute l0 ot l1 (negative log-likelihood)'
        end
    end

