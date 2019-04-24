 function [ c ] = LN( n )
        c0 = 2.865064;
        c = log2(c0);
        logTerm = log2(n);
        while logTerm > 0
            c = c + logTerm;
            logTerm = log2(logTerm);
        end
    end

