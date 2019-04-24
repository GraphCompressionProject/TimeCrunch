import config;
import mdl_base;
import mdl_structs;
import mdl_error;

from math import log,factorial;
from error import Error;
from graph import Graph;
from model import Model;
from temporal_model import *;
from copy import deepcopy;
from heapq import heappush, heappop, heapify
from collections import defaultdict

from mdl_base import *;
from mdl_structs import *;
from mdl_error import *;


### Our Encoding Starts Here ###

### Total Encoded Size
def L(G, M, errorEnc): 
    E = Error(G); # initially, everything is error, nothing is covered
    error_cost = 0;
 
   

    model_cost = LN(M.numStructs+1);    # encode number of structures we're encoding with
    model_cost += LwC(M.numStructs, M.numStrucTypes);            # encode the number per structure

    # encode the structure-type identifier per type
    if M.numFullCliques > 0 :
        model_cost += M.numFullCliques * log(M.numFullCliques / float(M.numStructs), 2);
    if M.numNearCliques  > 0 :
        model_cost += M.numNearCliques * log(M.numNearCliques / float(M.numStructs), 2);
    if M.numChains > 0 :
        model_cost += M.numChains * log(M.numChains / float(M.numStructs), 2);
    if M.numStars > 0 :
        model_cost += M.numStars * log(M.numStars / float(M.numStructs), 2);
    # off-diagonals
    if M.numFullOffDiagonals > 0 :
        model_cost += M.numFullOffDiagonals * log(M.numFullOffDiagonals / float(M.numStructs), 2);
    if M.numNearOffDiagonals > 0 :
        model_cost += M.numNearOffDiagonals * log(M.numNearOffDiagonals / float(M.numStructs), 2);
    # bipartite-cores
    if M.numBiPartiteCores > 0 :
        model_cost += M.numBiPartiteCores * log(M.numBiPartiteCores / float(M.numStructs), 2);
    if M.numNearBiPartiteCores > 0 :
        model_cost += M.numNearBiPartiteCores * log(M.numNearBiPartiteCores / float(M.numStructs), 2);
    if M.numJellyFishes > 0 :
        model_cost += M.numJellyFishes * log(M.numJellyFishes / float(M.numStructs), 2);
    if M.numCorePeripheries > 0 :
        model_cost += M.numCorePeripheries * log(M.numCorePeripheries / float(M.numStructs), 2);

    # encode the structures
    for struc in M.structs :
        if struc.isFullClique() :
            model_cost += LfullClique(struc,M,G,E);
        elif struc.isNearClique() :
            model_cost += LnearClique(struc,M,G,E);
        elif struc.isChain() :
            model_cost += Lchain(struc,M,G,E);
        elif struc.isStar() :
            model_cost += Lstar(struc,M,G,E);
        elif struc.isCorePeriphery() :
            model_cost += LcorePeriphery(struc,M,G,E);
        elif struc.isJellyFish() :
            model_cost += LjellyFish(struc,M,G,E);
        elif struc.isBiPartiteCore() :
            model_cost += LbiPartiteCore(struc,M,G,E);
        elif struc.isNearBiPartiteCore() :
            model_cost += LnearBiPartiteCore(struc,M,G,E);
        elif struc.isFullOffDiagonal() :
            model_cost += LfullOffDiagonal(struc,M,G,E);
        elif struc.isNearOffDiagonal() :
            model_cost += LnearOffDiagonal(struc,M,G,E);
    
    # encode the error
    error_cost += 0 if E.numCellsCovered == 0 else log(E.numCellsCovered, 2);    # encode number of additive Errors
    if ((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered > 0 :
        error_cost += log(((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered, 2);    # encode number of Errors
        
    if errorEnc == "NP" :
        error_cost += LErrorNaivePrefix(G,M,E);
    elif errorEnc == "NB" :
        error_cost += LErrorNaiveBinom(G,M,E);
    elif errorEnc == "TP" :
        error_cost += LErrorTypedPrefix(G,M,E);
    elif errorEnc == "TB" :
        error_cost += LErrorTypedBinom(G,M,E);
    
    total_cost = model_cost + error_cost;
    
    return (total_cost, model_cost, error_cost, E);

    
    
### Total Encoded Size for the greedy heuristic -- incrementally update the MDL cost
## for the newly added stucture 'struc'
def Lgreedy(G, M, errorEnc, time, struc, totalCostOld, Eold, model_cost_struct): 
    
    if time == 1:
        E = Error(G); # initially, everything is error, nothing is covered
        E.saveOld();
        # the cost for encoding each structure (to avoid recomputing it for the greedy updates)
        model_cost2 = 0;
    else :
        E = Eold;
        # the cost for encoding each structure separately 
        # Just update the up-to-now cost by adding the cost of the new structure
        model_cost2 = model_cost_struct;

    error_cost = 0;
    
    model_cost = LN(M.numStructs+1);    # encode number of structures we're encoding with
    model_cost += LwC(M.numStructs, M.numStrucTypes);            # encode the number per structure

    # encode the structure-type identifier per type
    if M.numFullCliques > 0 :
        model_cost += M.numFullCliques * log(M.numFullCliques / float(M.numStructs), 2);
    if M.numNearCliques  > 0 :
        model_cost += M.numNearCliques * log(M.numNearCliques / float(M.numStructs), 2);
    if M.numChains > 0 :
        model_cost += M.numChains * log(M.numChains / float(M.numStructs), 2);
    if M.numStars > 0 :
        model_cost += M.numStars * log(M.numStars / float(M.numStructs), 2);
    # off-diagonals
    if M.numFullOffDiagonals > 0 :
        model_cost += M.numFullOffDiagonals * log(M.numFullOffDiagonals / float(M.numStructs), 2);
    if M.numNearOffDiagonals > 0 :
        model_cost += M.numNearOffDiagonals * log(M.numNearOffDiagonals / float(M.numStructs), 2);
    # bipartite-cores
    if M.numBiPartiteCores > 0 :
        model_cost += M.numBiPartiteCores * log(M.numBiPartiteCores / float(M.numStructs), 2);
    if M.numNearBiPartiteCores > 0 :
        model_cost += M.numNearBiPartiteCores * log(M.numNearBiPartiteCores / float(M.numStructs), 2);
    if M.numJellyFishes > 0 :
        model_cost += M.numJellyFishes * log(M.numJellyFishes / float(M.numStructs), 2);
    if M.numCorePeripheries > 0 :
        model_cost += M.numCorePeripheries * log(M.numCorePeripheries / float(M.numStructs), 2);

    # encode the structures
    if struc.isFullClique() :
        model_cost2 += LfullClique(struc,M,G,E);
    elif struc.isNearClique() :
        model_cost2 += LnearClique(struc,M,G,E);
    elif struc.isChain() :
        model_cost2 += Lchain(struc,M,G,E);
    elif struc.isStar() :
        model_cost2 += Lstar(struc,M,G,E);
    elif struc.isCorePeriphery() :
        model_cost2 += LcorePeriphery(struc,M,G,E);
    elif struc.isJellyFish() :
        model_cost2 += LjellyFish(struc,M,G,E);
    elif struc.isBiPartiteCore() :
        model_cost2 += LbiPartiteCore(struc,M,G,E);
    elif struc.isNearBiPartiteCore() :
        model_cost2 += LnearBiPartiteCore(struc,M,G,E);
    elif struc.isFullOffDiagonal() :
        model_cost2 += LfullOffDiagonal(struc,M,G,E);
    elif struc.isNearOffDiagonal() :
        model_cost2 += LnearOffDiagonal(struc,M,G,E);
    
    # encode the error
    error_cost += 0 if E.numCellsCovered == 0 else log(E.numCellsCovered, 2);    # encode number of additive Errors
    if ((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered > 0 :
        error_cost += log(((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered, 2);    # encode number of Errors
 
    if errorEnc == "NP" :
        error_cost += LErrorNaivePrefix(G,M,E);
    elif errorEnc == "NB" :
        error_cost += LErrorNaiveBinom(G,M,E);
    elif errorEnc == "TP" :
        error_cost += LErrorTypedPrefix(G,M,E);
    elif errorEnc == "TB" :
        error_cost += LErrorTypedBinom(G,M,E);
    
    total_cost = model_cost + model_cost2 + error_cost;
    model_cost_total = model_cost + model_cost2;    

    return (total_cost, model_cost_total, model_cost2, error_cost, E);

### Total Encoded Size for temporal model
def Ltemporal(G_arr, TM, errorEnc):
    
    #first, determine total cost of storing the model itself.  after doing this, compute error costs for each timestep
    
    E = Error(G_arr[0]); # just use any error object for now - we don't care about error right now - only considering the model cost
    
    #these are counted in the following line, so we comment them out
    #temporal_model_cost_0 = LN(TM.numStructs+1);    # encode number of structures we're encoding with
    #temporal_model_cost_0 += LwC(TM.numStructs, TM.numStrucTypes);            # encode the number per structure

    temporal_model_cost_0  = L(G_arr[0], TM, errorEnc)[1]; #cost of encoding just the structures (model) without time information
    ntimesteps = len(G_arr);
    for struc in TM.structs: #for each struct, encode temporal information
        temporal_type = struc.getTemporalType();
	time_costs = get_tsig(struc.getTimesteps(), ntimesteps);
        if temporal_type == "o":
            temporal_model_cost_0 += time_costs[1];
        elif temporal_type == "c":
            temporal_model_cost_0 += time_costs[5];
        elif temporal_type == "p":
            temporal_model_cost_0 += time_costs[3]; #encode start, step size, and how many steps to take
        elif temporal_type == "r":
            temporal_model_cost_0 += time_costs[4]; #encode start and how many steps to take (step size is always 1)
        elif temporal_type == "f":
	    temporal_model_cost_0 += time_costs[2];	
    # encode the structure-type identifier per type
    #if TM.numFullCliques > 0 :
    #    temporal_model_cost_0 += TM.numFullCliques * log(TM.numFullCliques / float(TM.numStructs), 2);
    #if TM.numNearCliques  > 0 :
    #    temporal_model_cost_0 += TM.numNearCliques * log(TM.numNearCliques / float(TM.numStructs), 2);
    #if TM.numChains > 0 :
    #    temporal_model_cost_0 += TM.numChains * log(TM.numChains / float(TM.numStructs), 2);
    #if TM.numStars > 0 :
    #    temporal_model_cost_0 += TM.numStars * log(TM.numStars / float(TM.numStructs), 2);
    ### off-diagonals
    #if TM.numFullOffDiagonals > 0 :
    #    temporal_model_cost_0 += TM.numFullOffDiagonals * log(TM.numFullOffDiagonals / float(TM.numStructs), 2);
    #if TM.numNearOffDiagonals > 0 :
    #    temporal_model_cost_0 += TM.numNearOffDiagonals * log(TM.numNearOffDiagonals / float(TM.numStructs), 2);
    ## bipartite-cores
    #if TM.numBiPartiteCores > 0 :
    #    temporal_model_cost_0 += TM.numBiPartiteCores * log(TM.numBiPartiteCores / float(TM.numStructs), 2);
    #if TM.numNearBiPartiteCores > 0 :
    #    temporal_model_cost_0 += TM.numNearBiPartiteCores * log(TM.numNearBiPartiteCores / float(TM.numStructs), 2);
    #if TM.numJellyFishes > 0 :
    #    temporal_model_cost_0 += TM.numJellyFishes * log(TM.numJellyFishes / float(TM.numStructs), 2);
    #if TM.numCorePeripheries > 0 :
    #    temporal_model_cost_0 += TM.numCorePeripheries * log(TM.numCorePeripheries / float(TM.numStructs), 2);
    #all model costs have been taken care of now.  we now want to compute error costs for each timestep
    total_error_cost = 0;
    E_arr = [];
    error_cost_arr = [];
    model_cost_arr = [];
    for tstep in range(0,len(G_arr)):
        G = G_arr[tstep];
        model_cost = 0;
        error_cost = 0;
        E = Error(G); #initially, everything is error

        #initialize a new "relevant" model which only includes the structures which appear at this timestep, since we don't need to encode all temporal structures for each slice
        RTM = TModel();
        #for encoding structures within a single timestep, we just need non-temporal structure types (just fc, nc, ch, st, bc and nb) so we divide by the temporal types
        RTM.numStrucTypes = RTM.numStrucTypes/5;
        relv_structs = [TM.structs[k] for k in range(0,TM.numStructs) if tstep+1 in TM.structs[k].getTimesteps()]
        #populate this new model with the appropriate structures
        for struc in relv_structs:
            RTM.addStructure(struc);

        # encode the structures
        for struc in RTM.structs :
            if struc.isFullClique() :
                model_cost += LfullClique(struc,RTM,G,E);
            elif struc.isNearClique() :
                model_cost += LnearClique(struc,RTM,G,E);
            elif struc.isChain() :
                model_cost += Lchain(struc,RTM,G,E);
            elif struc.isStar() :
                model_cost += Lstar(struc,RTM,G,E);
            elif struc.isCorePeriphery() :
                model_cost += LcorePeriphery(struc,RTM,G,E);
            elif struc.isJellyFish() :
                model_cost += LjellyFish(struc,RTM,G,E);
            elif struc.isBiPartiteCore() :
                model_cost += LbiPartiteCore(struc,RTM,G,E);
            elif struc.isNearBiPartiteCore() :
                model_cost += LnearBiPartiteCore(struc,RTM,G,E);
            elif struc.isFullOffDiagonal() :
                model_cost += LfullOffDiagonal(struc,RTM,G,E);
            elif struc.isNearOffDiagonal() :
                model_cost += LnearOffDiagonal(struc,RTM,G,E);
	
        # encode the error
        error_cost += 0 if E.numCellsCovered == 0 else log(E.numCellsCovered, 2);    # encode number of additive Errors from the relevant model
        if ((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered > 0 :
            error_cost += log(((G.numNodes * G.numNodes - G.numNodes) / 2) - E.numCellsCovered, 2);    # encode number of Errors from the relevant model

        if errorEnc == "NP" :
            error_cost += LErrorNaivePrefix(G,RTM,E);
        elif errorEnc == "NB" :
            error_cost += LErrorNaiveBinom(G,RTM,E);
        elif errorEnc == "TP" :
            error_cost += LErrorTypedPrefix(G,RTM,E);
        elif errorEnc == "TB" :
            error_cost += LErrorTypedBinom(G,RTM,E);

        E_arr.append(E);
        error_cost_arr.append(error_cost);
        model_cost_arr.append(model_cost);
        total_error_cost += error_cost;
    total_cost = temporal_model_cost_0 + total_error_cost;
    return (total_cost, temporal_model_cost_0, total_error_cost, model_cost_arr, error_cost_arr, E_arr);


def Ltemporal_greedy(G_arr, TM, errorEnc):
	
	#start with empty model
	CM = TModel();
	(total_cost_orig, model_cost_orig, error_cost_orig, model_cost_arr_orig, error_cost_arr_orig, E_arr) = Ltemporal(G_arr, CM, errorEnc);
	kept_strucs = [];
	prev_total_cost = total_cost_orig; #total_cost_orig;
	prev_model_cost = 0; #model_cost_orig;

	print "initial total cost ", prev_total_cost
	print "---------------------------------------------"
	
	#now that we have the total with no structures, start adding structures to see if the total increases or decreases
	
	for cid in range(0,TM.numStructs):
		
		candidate = TM.structs[cid]

		#initialize the model cost of storing the structures up until now (0 at the start)
		model_cost = prev_model_cost;
		print "model cost initialized to ",model_cost
		
		#save all the errors that we start with in case the structure isn't worth encoding 
		for E in E_arr:
			E.saveOld();

		#add candidate
		CM.addStructure(candidate);
		print "trying out structure ", cid
		#encode how many structures and how many of each
		model_cost += LN(CM.numStructs+1);
		model_cost += LwC(CM.numStructs, CM.numStrucTypes);
		# encode the structure-type identifier per type
	
    		if CM.numFullCliques > 0 :
        		model_cost += CM.numFullCliques * log(CM.numFullCliques / float(CM.numStructs), 2);
		if CM.numNearCliques  > 0 :
		        model_cost += CM.numNearCliques * log(CM.numNearCliques / float(CM.numStructs), 2);
		if CM.numChains > 0 :
		        model_cost += CM.numChains * log(CM.numChains / float(CM.numStructs), 2);
		if CM.numStars > 0 :
		        model_cost += CM.numStars * log(CM.numStars / float(CM.numStructs), 2);
		# off-diagonals
		if CM.numFullOffDiagonals > 0 :
		        model_cost += CM.numFullOffDiagonals * log(CM.numFullOffDiagonals / float(CM.numStructs), 2);
		if CM.numNearOffDiagonals > 0 :
		        model_cost += CM.numNearOffDiagonals * log(CM.numNearOffDiagonals / float(CM.numStructs), 2);
		# bipartite-cores
		if CM.numBiPartiteCores > 0 :
		        model_cost += CM.numBiPartiteCores * log(CM.numBiPartiteCores / float(CM.numStructs), 2);
		if CM.numNearBiPartiteCores > 0 :
		        model_cost += CM.numNearBiPartiteCores * log(CM.numNearBiPartiteCores / float(CM.numStructs), 2);
		if CM.numJellyFishes > 0 :
		        model_cost += CM.numJellyFishes * log(CM.numJellyFishes / float(CM.numStructs), 2);
		if CM.numCorePeripheries > 0 :
		        model_cost += CM.numCorePeripheries * log(CM.numCorePeripheries / float(CM.numStructs), 2);
		ntimesteps = len(G_arr);

		#used to compute the temporal + model cost for this particular candidate structure
		model_cost_candidate = 0;
		#encode the structure's temporal cost
    		temporal_type = candidate.getTemporalType();
        	time_costs = get_tsig(candidate.getTimesteps(), ntimesteps);
        	if temporal_type == "o":
			model_cost_candidate += time_costs[1];
        	elif temporal_type == "c":
            		model_cost_candidate += time_costs[5];
        	elif temporal_type == "p":
            		model_cost_candidate += time_costs[3]; #encode start, step size, and how many steps to take
        	elif temporal_type == "r":
            		model_cost_candidate += time_costs[4]; #encode start and how many steps to take (step size is always 1)
        	elif temporal_type == "f":
            		model_cost_candidate += time_costs[2];
		
		# encode the structure's model cost in terms of edges
        	if candidate.isFullClique() :
			tstep = candidate.getTimesteps()[0];
			d = model_cost_candidate;
			model_cost_candidate += LfullClique(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
			print model_cost_candidate - d	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LfullClique(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isNearClique() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LnearClique(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LnearClique(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isChain() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += Lchain(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				Lchain(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isStar() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += Lstar(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
			other_tsteps = candidate.getTimesteps()[1:];
			for tstep in other_tsteps:
				Lstar(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isCorePeriphery() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LcorePeriphery(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LcorePeriphery(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isJellyFish() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LjellyFish(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LjellyFish(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isBiPartiteCore() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LbiPartiteCore(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LbiPartiteCore(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isNearBiPartiteCore() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LnearBiPartiteCore(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LnearBiPartiteCore(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isFullOffDiagonal() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LfullOffDiagonal(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LfullOffDiagonal(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
            	elif candidate.isNearOffDiagonal() :
                	tstep = candidate.getTimesteps()[0];
			model_cost_candidate += LnearOffDiagonal(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
	
			other_tsteps = candidate.getTimesteps()[1:];

			for tstep in other_tsteps:
				LnearOffDiagonal(candidate, CM, G_arr[tstep-1], E_arr[tstep-1]);
			
		#now we need to compute error costs for each graph
		error_cost = 0;
		for i in range(0,ntimesteps):
			error_cost += 0 if E_arr[i].numCellsCovered == 0 else log(E_arr[i].numCellsCovered, 2);    # encode number of additive Errors from the relevant model
			if ((G_arr[i].numNodes * G_arr[i].numNodes - G_arr[i].numNodes) / 2) - E_arr[i].numCellsCovered > 0 :
				error_cost += log(((G_arr[i].numNodes * G_arr[i].numNodes - G_arr[i].numNodes) / 2) - E_arr[i].numCellsCovered, 2);    # encode number of Errors from the relevant model

		        if errorEnc == "NP" :
				error_cost += LErrorNaivePrefix(G_arr[i], CM, E_arr[i]);
			elif errorEnc == "NB" :
				error_cost += LErrorNaiveBinom(G_arr[i], CM, E_arr[i]);
			elif errorEnc == "TP" :
				error_cost += LErrorTypedPrefix(G_arr[i], CM, E_arr[i]);
			elif errorEnc == "TB" :
				error_cost += LErrorTypedBinom(G_arr[i], CM, E_arr[i]);

		#compute the total cost by adding the parts of the other parts of the model to the structure costs (time and model)
		total_cost = model_cost + model_cost_candidate + error_cost;
		
		#if the new cost is greater than the old one, discard structure and recover errors.  else set prev_total_cost to total_cost and continue with future structures
		if total_cost > prev_total_cost:
			print "previous total cost ", prev_total_cost
			print "new total cost ", total_cost
			print "discarding structure..."
			print "---------------------------------------------"
			#discard structure
			CM.rmStructure(candidate);
			#for each error object, recover the old state of errors, before we had added this structure
			for E in E_arr:
				E.recoverOld();
		else:
			kept_strucs.append(cid);
			print "previous total cost ", prev_total_cost
			print "new total cost ", total_cost
			print "kept the structure..."
			print "---------------------------------------------"
			prev_total_cost = total_cost; #we have a new cost to improve upon in the next iteration
			prev_model_cost += model_cost_candidate; #add the previous cost to the candidate cost thus far.  this will build up across iterations to contain only structure costs (time and model)
	return kept_strucs;	

#auxiliary functions for determining temporal signature (diff, mean, stdev)
def diff(a):
        #compute first order differences
        return [a[i+1]-a[i] for i in range(0,len(a)-1)];

def mean(a):
        #compute arithmetic mean
        return float(sum(a))/len(a);

def stdev(a):
        #compute standard deviation
        avg_a = mean(a);
        var_comp = map(lambda x: (x - avg_a)**2, a);
        return math.sqrt(mean(var_comp));
    
def get_tsig(tsteps, ntot):
        scores = ['?'];
        scores.extend([float("inf")] * 5); # [oneshot, flickering, periodic, range, constant] 

        if len(tsteps) == 1:
                scores[1] = log(ntot, 2);
                scores[0] = 'o';
                return scores;
        elif len(tsteps) == ntot:
                scores[5] = 0;
                scores[0] = 'c';
                return scores;

        d_tsteps = diff(tsteps);
        pdcity = mean(d_tsteps);

        #compute cost of periodic/range (special case) encoding 
        projected = [round(tsteps[0]+(i*pdcity)) for i in range(0,len(tsteps))]; #need an equal length series to compare two 'models'

        model_cost_periodic = log(choose(ntot, 2), 2) + LN(len(tsteps));
        errors = [(tsteps[i] - projected[i]) for i in range(0,len(projected))];
	#count unique elements in errors to encode w/ prefix codes	

	total = len(errors);
	hist = defaultdict(int);
	for i in errors:
		hist[i] += 1;
	
	vals = hist.keys();
	error_cost_periodic = 0;
	if len(vals) >= 2:
		for i in vals:
			error_cost_periodic += (-hist[i] * log(float(hist[i])/total, 2)); 
	
        total_cost_periodic = model_cost_periodic + error_cost_periodic;

        total_cost_flickering = log(choose(ntot, len(tsteps)), 2) + LN(len(tsteps));
        if total_cost_flickering < total_cost_periodic:
                scores[2] = total_cost_flickering;
                scores[0] = 'f';
        else:
                if pdcity >= 1.5:
                        scores[3] = total_cost_periodic;
                        scores[0] = 'p'
                else:
                        scores[4] = total_cost_periodic;
                        scores[0] = 'r'
        return scores;
    
