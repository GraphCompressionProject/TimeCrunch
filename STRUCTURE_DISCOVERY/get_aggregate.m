function [model_struc] = get_aggregate(filelist, minSize)


addpath('STRUCTURE_DISCOVERY');
addpath('STRUCTURE_DISCOVERY/VariablePrecisionIntegers/VariablePrecisionIntegers');
addpath('gaimc/');

inpfile = fopen(filelist);
fnames = textscan(inpfile,'%s');
ntimesteps = numel(fnames{1});

aggregate_graph = [];

for i = 1:ntimesteps
	input_file = fnames{1}{i};
	disp(input_file);
	temp = load(input_file);
	blah = find(temp == 0);
	if(numel(blah) > 0)
		quit;
	end
 	if(size(temp,2) < 3)
 	       temp(:,3) = 1;
 	end
	aggregate_graph = [aggregate_graph; temp];
end
size(aggregate_graph)
min(min(aggregate_graph))
max(max(aggregate_graph))	
orig = spconvert(aggregate_graph);
if(size(orig,1) ~= size(orig,2))
	orig(max(size(orig)),max(size(orig))) = 0;
end

n = max(size(orig));

orig_sym = orig + orig';
orig_sym = orig_sym~=0;
orig_sym_nodiag = orig_sym - diag(diag(orig_sym));
[i,j,~] = find(orig_sym_nodiag);
dlmwrite('aggregateGraph.txt',[i,j],'precision','%d','delimiter',' ');


global model;
model = struct('code', {}, 'edges', {}, 'nodes1', {}, 'nodes2', {}, 'benefit', {}, 'benefit_notEnc', {});
global model_idx;
model_idx = 0;
model_struc = SlashBurnEncode( orig_sym_nodiag, 4, '.', false, false, minSize, 'aggregate');

end	

