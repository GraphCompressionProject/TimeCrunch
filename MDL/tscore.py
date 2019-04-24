#!/usr/local/bin/python2.6

import sys
import os
import config
import math

from pprint import pprint

from time import time

from mdl import *;
from error import Error;
from graph import Graph;
from model import *;

from temporal_model import *;

#from description_length import *;


if len(sys.argv) <= 1 :
    print 'at least: <graph-filelist> [model.model] [-pC] [-lC] [-pE] [-lE] [-e{NP,NB,TP,TB}]';
    print ' optional argument model.model = file to read model from, otherwise only empty model';
    print ' optional argument -vX    = verbosity (1, 2, or 3)';
    print ' optional argument -pG    = plot Graph adjacency matrix';
    print ' optional argument -pC    = plot Cover matrix';
    print ' optional argument -pE    = plot Error matrix';
    print ' optional argument -lC    = list Cover entries';
    print ' optional argument -lE    = list Error entries';
    print ' optional argument -eXX   = encode error resp. untyped using prefix (NP), or';
    print '                            binomial (NB) codes, or using typed';
    print '                            prefix (TP) or binomial (TB, default) codes';
    exit();

if (len(sys.argv) > 1 and ("-v1" in sys.argv)) :
    config.optVerbosity = 1;
elif (len(sys.argv) > 1 and ("-v2" in sys.argv)) :
    config.optVerbosity = 2;
if (len(sys.argv) > 1 and ("-v3" in sys.argv)) :
    config.optVerbosity = 3;

t0 = time()

gFilename_list = sys.argv[1];
graphname_list = [x.rstrip() for x in open(gFilename_list,'r').readlines()]
graph_list = list()

nnodes = 0;

for graph_name in graphname_list:
	g = Graph();
	g.load(graph_name);
	graph_list.append(g);
	nnodes = max(nnodes,g.numNodes);

print "nnodes: " + str(nnodes)

for tstep in range(0,len(graph_list)):
	print graphname_list[tstep];
	g = graph_list[tstep];
	g.numNodes = nnodes;
	currentEdgeLen = len(g.edges);
	g.edges.extend([frozenset() for x in range(currentEdgeLen,nnodes)]);
		


temp_model_file = sys.argv[2];

#instantiate empty temporal model
TM = TModel();

#set error type
errorEnc = config.optDefaultError;

print "   \t" + "L(G,M)" + "\tL(M)" + "\tL(E)";

#get costs without using a model
(total_cost, temporal_model_cost, total_error_cost, model_cost_arr, error_cost_arr, E_arr) = Ltemporal(graph_list, TM, errorEnc);
temporal_model_cost += len(graph_list) * log(len(graph_list), 2);
total_cost += len(graph_list) * log(len(graph_list), 2);
print "M_0:\t" + '%.0f' % total_cost + "\t" + '%.0f' % temporal_model_cost + "\t" + '%.0f' % total_error_cost + "\t";

#load the model
TM.load(temp_model_file);

#get costs, this time with a model
(total_cost, temporal_model_cost, total_error_cost, model_cost_arr, error_cost_arr, E_arr) = Ltemporal(graph_list, TM, errorEnc);

print "M_x:\t" + '%0.f' % total_cost + "\t" + '%0.f' % temporal_model_cost + "\t" + '%0.f' % total_error_cost + "\t";
