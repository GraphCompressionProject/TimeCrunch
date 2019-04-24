 function [ c_err ] = Lnu_opt( E )
        % E has two entries: # of included edges, # of excluded edges
        Einc = E(1);
        Eexc = E(2);
        c_err = LN( Einc ) + ...
            Einc * NLL( Einc, Eexc, 1) + ...
            Eexc * NLL( Einc, Eexc, 0);
    end
