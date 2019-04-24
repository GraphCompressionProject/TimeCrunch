model_filelist = 'path/to/file/containing/ordered/list/of/model/files';
graph_filelist = 'path/to/file/containing/ordered/list/of/graph/files';
nstrucs_wanted = 5000;
thresh = 0.5;
outfile = ['timecrunch-',num2str(thresh),'.tmodel'];

stars_edges = [];
bips_edges = [];
cliqs_edges = [];
chains_edges = [];

star_hubs = [];

star_ct = 0;
bip_ct = 0;
cliq_ct = 0;
chain_ct = 0;

star_times = [];
bip_times = [];
cliq_times = [];
chain_times = [];

tstep = 0;

star_orig = cell(1);
bip_orig = cell(1);
cliq_orig = cell(1);
chain_orig = cell(1);

flid = fopen(model_filelist);
mfn = fgetl(flid);

while ischar(mfn)
	disp(['reading ',mfn,'...']);
	fid = fopen(mfn);
	tstep = tstep + 1;

	tline = fgetl(fid);
	while ischar(tline)
		pieces = strread(tline,'%s','delimiter',', ')';
		str_nodes = pieces(2:end);
		edges = cellfun(@str2num,str_nodes);
		
		if pieces{1} == 'st'
			star_ct = star_ct + 1;
			stars_edges = [stars_edges; horzcat(ones(numel(edges), 1)*star_ct, edges',ones(numel(edges), 1))];
			star_times = [star_times; tstep];
			star_hubs = [star_hubs; edges(1)];
			star_orig{star_ct} = tline;
		elseif pieces{1} == 'fc' | pieces{1} == 'nc'
			cliq_ct = cliq_ct + 1;
			cliqs_edges = [cliqs_edges; horzcat(ones(numel(edges), 1)*cliq_ct, edges',ones(numel(edges), 1))];
			cliq_times = [cliq_times; tstep];
			cliq_orig{cliq_ct} = tline;
		elseif pieces{1} == 'bc' | pieces{1} == 'nb'
			bip_ct = bip_ct + 1;
			bips_edges = [bips_edges; horzcat(ones(numel(edges), 1)*bip_ct, edges',ones(numel(edges), 1))];
			bip_times = [bip_times; tstep];
			bip_orig{bip_ct} = tline;
		elseif pieces{1} == 'ch'
			chain_ct = chain_ct + 1;
			chains_edges = [chains_edges; horzcat(ones(numel(edges), 1)*chain_ct, edges',ones(numel(edges), 1))];
			chain_times = [chain_times; tstep];
			chain_orig{chain_ct} = tline;
		end
		tline = fgetl(fid);
	end
	fclose(fid);
	mfn = fgetl(flid);
end
fclose(flid);

disp('all files read');

orig_lines = {star_orig, cliq_orig, bip_orig, chain_orig};
status = [0,0,0,0];
if numel(stars_edges) > 0
	star_mat = spconvert(stars_edges);
	status(1) = 1;
else
	star_mat = sparse(1,1);
end

if numel(cliqs_edges) > 0
	cliq_mat = spconvert(cliqs_edges);
	status(2) = 1;
else
	cliq_mat = sparse(1,1);
end

if numel(bips_edges) > 0
	bip_mat = spconvert(bips_edges);
	status(3) = 1;
else
	bip_mat = sparse(1,1);
end

if numel(chains_edges) > 0
	chain_mat = spconvert(chains_edges);
	status(4) = 1;
else
	chain_mat = sparse(1,1);
end

%handle stars, cliques, bipartites, and chains
smms = {star_mat, cliq_mat, bip_mat, chain_mat};

disp('smms created');

ntstrucs = 0;
tstrucs = struct('tsig','','struc_str','','tsteps',[],'enc_benefit',0);

while ntstrucs <= nstrucs_wanted
	if nnz(status) == 0
		break
	end

	setidx = randsample(find(status), 1);

	mat = smms{setidx};
	wtime = 0;
	disp(['ntstrucs so far: ',num2str(ntstrucs)]);	
		
	consecutive_fails = 0;

	if setidx ~= 1

		%rank-1 SVD, and sort based on projection magnitude
		tic; 
		[u,s,v] = svds(mat,1);
		[~,srt] = sort(abs(u(:,1)), 'descend');

		%initialize matches vector, start with the base
		matches = [srt(1)];
		base = mat(srt(1),:);
		
		%for each other structure, try to match, and if it doesn't match 10 times in a row, move on to next iteration
		for j = 2:size(mat,1)
			candidate = mat(srt(j),:);

			if sum(base & candidate)/sum(base | candidate) < thresh
				consecutive_fails = consecutive_fails + 1;
				if consecutive_fails == 10
					break;
				end
				continue;
			end
			consecutive_fails = 0;
			matches = [matches; srt(j)];
		end

		%what to do with matches
		matches = sort(matches);
		
		if setidx == 2
			struc_tsteps = unique(cliq_times(matches))';
		elseif setidx == 3
			struc_tsteps = unique(bip_times(matches))';
		elseif setidx == 4
			struc_tsteps = unique(chain_times(matches))';
		end
	
		tsig = get_tsig(struc_tsteps, tstep);
	
		%clear matched structures so we don't reuse them
		for t = 1:numel(matches)
			mat(matches(t),:) = 0;
		end
	
		smms{setidx} = mat;
		
		%if we don't have any structures left to process, move on to the next structure set
		if nnz(mat) < 1
			status(setidx) = 0;
		end
	
		wtime = wtime + toc;

		line_set = orig_lines{setidx};
		base_line = line_set{srt(1)};
		
		ntstrucs = ntstrucs + 1;
		t_strucs(ntstrucs).tsig = tsig;
		t_strucs(ntstrucs).struc_str = base_line;
		t_strucs(ntstrucs).tsteps = struc_tsteps; 	
	else
		if nnz(star_hubs) > 1
			hub = star_hubs(randsample(find(star_hubs), 1));
		else
			hub = star_hubs(find(star_hubs));
		end

		star_ids = find(star_hubs == hub);
		mat2 = mat(star_ids, :);
        
		tic;
	        %rank-1 SVD, and sort based on projection magnitude 
                [u,s,v] = svds(mat2,1);
                [~,srt] = sort(abs(u(:,1)), 'descend');

                %initialize matches vector, start with the base
                matches = [srt(1)];
                base = mat2(srt(1),:);

                %for each other structure, try to match, and if it doesn't match 10 times in a row, move on to next iteration
                for j = 2:size(mat2,1)
                        candidate = mat2(srt(j),:);

                        if sum(base & candidate)/sum(base | candidate) < thresh
                                consecutive_fails = consecutive_fails + 1;
                                if consecutive_fails == 10
                                        break;
                                end
                                continue;
                        end
                        consecutive_fails = 0;
                        matches = [matches; srt(j)];
                end

		matches = star_ids(matches);	
		matches = sort(matches);
		struc_tsteps = unique(star_times(matches))';
		
		tsig = get_tsig(struc_tsteps, tstep);

		for t = 1:numel(matches)
			mat(matches(t),:) = 0;
			star_hubs(matches(t)) = 0;	
		end	

		smms{setidx} = mat;
		
		if nnz(mat) < 1
			status(setidx) = 0;
		end

		wtime = wtime + toc;

		line_set = orig_lines{setidx};
                base_line = line_set{star_ids(srt(1))};

                ntstrucs = ntstrucs + 1;
                t_strucs(ntstrucs).tsig = tsig;
                t_strucs(ntstrucs).struc_str = base_line;
                t_strucs(ntstrucs).tsteps = struc_tsteps;
	
	
	end
		
end

%now we need to rank the structures
graphs = cell(1,tstep);
gflid = fopen(graph_filelist);

ctr = 0;
nnodes = 0;

gfn = fgetl(gflid);

%read all the graphs, and convert to undirected, sparse matrices
while ischar(gfn)
	disp(['reading ',gfn,'...']);
	g = load(gfn);
	%load each sparse graph	
	if numel(g) > 0
		g(:,3) = 1;
		g = spconvert(g);
	else
		g = sparse(1,1);
	end
	
	%binarize and remove diagonal
	g(max(size(g)), max(size(g))) = 0;
	g = g + g';
	g = g - diag(diag(g));
	g(g~=0) = 1;
	
	%keep track of max nodecount
	nnodes = max(nnodes, max(size(g)));
	
	ctr = ctr + 1;
	graphs{ctr} = g;
	gfn = fgetl(gflid);
end

fclose(gflid);

%augment size of each graph to max actual size
for k = 1:ctr
	graphs{k}(nnodes,nnodes) = 0;
end

disp('all graphs read and processed');
%compute encoding benefit for each struct
tic;
for c = 1:ntstrucs
	base_str = t_strucs(c).struc_str;
	pieces = strsplit(base_str,', ');
	struc_type = base_str(1:2);
	nodes = [];
	lnodes = [];
	rnodes = [];

	%string processing to get nodeset, and left/right node split if necessary
	if numel(pieces) == 1
		pieces = strread(pieces{1},'%s','delimiter',' ')';
		nodes = cellfun(@str2num,pieces(2:end));
	elseif numel(pieces) == 2
		pieces2 = strread(pieces{1},'%s','delimiter',' ')';
		if numel(pieces2) == 2
			lnodes = [str2num(pieces2{2})];
			rnodes = [str2num(pieces{2})];
			nodes = [lnodes rnodes];
		else
			lnodes = cellfun(@str2num, pieces2(2:end)); 
			rnodes = str2num(pieces{2});			
			nodes = [lnodes rnodes];
		end
	end
	
	n_sub = numel(nodes);
	nleft = numel(lnodes);
	nright = numel(rnodes);
	tsteps = t_strucs(c).tsteps;
	
	%encoding benefit for t. strucs is computed as sum(cost_as_error/cost_with_model)
	encoding_benefit = 0;

	%for each timestep, get the subgraph in the appropriate time graph and compute errors
	disp(['computing local encoding benefit for structure ',num2str(c),' of ',num2str(ntstrucs)])
	for u = 1:numel(tsteps)
		i = tsteps(u);
		Asmall = graphs{i}(nodes, nodes);

		%get cost as an error matrix
		E = [nnz(Asmall) n_sub^2 - n_sub - nnz(Asmall)];
		cost_aserror = compute_encodingCost('err', 0, 0, E);
		
		%get model+error cost based on which kind of structure this is
		if struc_type == 'st'
			E(1) = 2*(n_sub - 1 - nnz(graphs{i}(lnodes,rnodes))) + nnz(graphs{i}(rnodes, rnodes));
			E(2) = n_sub^2 - n_sub - E(1);
			MDLcost = compute_encodingCost('st', nnodes, n_sub, E);
		elseif struc_type == 'fc' 
			E(1) = nnz(Asmall);
			E(2) = n_sub^2 - n_sub - E(1);
			MDLcost = compute_encodingCost('fc', nnodes, n_sub, E);
		elseif struc_type == 'nc'
			MDLcost = compute_encodingCost('nc', nnodes, n_sub, Asmall);
		elseif struc_type == 'bc'
			E(1) = 2*(nleft*nright - nnz(graphs{i}(lnodes,rnodes))) + nnz(graphs{i}(lnodes,lnodes)) + nnz(graphs{i}(rnodes,rnodes));
			E(2) = n_sub^2 - n_sub - E(1);
			MDLcost = compute_encodingCost('bc', nnodes, nleft, E, nright); 
		elseif struc_type == 'nb'
			E(1) = nnz(graphs{i}(lnodes, lnodes)) + nnz(graphs{i}(rnodes, rnodes));
                        E(2) = n_sub^2 - n_sub - E(1);		
			nb_edges = [nnz(graphs{i}(lnodes,rnodes)), nleft*nright-nnz(graphs{i}(lnodes,rnodes))];
			MDLcost = compute_encodingCost('nb', nnodes, nleft, E, nright, nb_edges);
		elseif struc_type == 'ch'
			missing = 0;
			existing = 0;
			for j = 1:n_sub-1
				if Asmall(j,j+1) == 0
					missing = missing+1;
				else
					existing = existing+1;
				end
			end
			E(1) = 2*missing + (nnz(Asmall) - 2*existing);
			E(2) = n_sub^2 - n_sub - E(1);
			MDLcost = compute_encodingCost('ch', nnodes, n_sub, E);
		end

		encoding_benefit = encoding_benefit + (cost_aserror/MDLcost);
	end 

	t_strucs(c).enc_benefit = encoding_benefit;
end

[~,order] = sort([t_strucs(:).enc_benefit], 'descend');
t_strucs_ordered = t_strucs(order);

s_time = toc;

fsout = fopen(outfile,'w');

for c = 1:numel(t_strucs_ordered)
	tm_str = [t_strucs_ordered(c).tsig{1}, t_strucs_ordered(c).struc_str];
	fprintf(fsout, '%s, ',tm_str);
	if numel(t_strucs_ordered(c).tsteps) > 1
		fprintf(fsout,'%d ',t_strucs_ordered(c).tsteps(1:end-1));
	end
	fprintf(fsout,'%d\n',t_strucs_ordered(c).tsteps(end));
end

fclose(fsout);

disp(['computation time ',num2str(wtime)])
disp(['evaluating and sorting time ',num2str(s_time)]);

