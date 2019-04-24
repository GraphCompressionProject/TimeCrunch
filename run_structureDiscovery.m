function [] = run_structureDiscovery(input_file) 

 unweighted_graph = input_file;
 output_model_greedy = 'path/to/output/models'

 addpath('STRUCTURE_DISCOVERY');
 addpath('STRUCTURE_DISCOVERY/VariablePrecisionIntegers/VariablePrecisionIntegers');
 addpath('gaimc/');
 orig = load(input_file);
 if(size(orig,2) < 3)
	orig(:,3) = 1;
 end
 orig = spconvert(orig);
 orig(max(size(orig)),max(size(orig))) = 0;
 orig_sym = orig + orig';
 orig_sym = orig_sym~=0;
 orig_sym_nodiag = orig_sym - diag(diag(orig_sym));
 disp('file LOADED')
 disp('==== Running VoG for structure discovery ====')
 global model; 
 model = struct('code', {}, 'edges', {}, 'nodes1', {}, 'nodes2', {}, 'benefit', {}, 'benefit_notEnc', {});
 global model_idx;
 model_idx = 0;
 SlashBurnEncode( orig_sym_nodiag, 2, output_model_greedy, false, false, 5, unweighted_graph);
end

% quit
