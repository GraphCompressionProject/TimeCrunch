function [ costs ] = get_tsig(tsteps, ntot)
	costs = {'?',Inf,Inf,Inf,Inf,Inf};
	if numel(tsteps) == 1
		costs{2} = log2(ntot);
		costs{1} = 'o';
		return;
	elseif numel(tsteps) == ntot
		costs{6} = 0;
		costs{1} = 'c';
		return;
	end

	d_tsteps = diff(tsteps);
	projected = round(linspace(tsteps(1), tsteps(end), numel(tsteps)));
	pdcity = projected(2) - projected(1);

	model_cost_periodic = l2cnk(ntot, 2) + LN(numel(tsteps));
	errors = tsteps - projected;
	[counts, vals] = hist(errors,unique(errors));
	probs = counts / sum(counts);

	error_cost_periodic = 0;	
	if numel(vals) < 2
		error_cost_periodic = 0;
	else
		for i = 1:numel(vals)
			count = counts(i);
			prob = probs(i);
			error_cost_periodic = error_cost_periodic + (-count * log2(prob));
		end		
	end
	total_cost_periodic = model_cost_periodic + error_cost_periodic;
	total_cost_flickering = l2cnk(ntot, numel(tsteps)) + LN(numel(tsteps));
	if total_cost_flickering < total_cost_periodic
		costs{3} = total_cost_flickering;
		costs{1} = 'f';
	else
		if pdcity >= 1.5
			costs{4} = total_cost_periodic;
			costs{1} = 'p';
		else
			costs{5} = total_cost_periodic;
			costs{1} = 'r';
		end
	end
end

 function [ c ] = LN( n )
        c0 = 2.865064;
        c = log2(c0);
        logTerm = log2(n);
        while logTerm > 0
            c = c + logTerm;
            logTerm = log2(logTerm);
        end
    end

