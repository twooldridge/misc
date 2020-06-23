#!/usr/bin/env python -u
import re, sys, os, getopt
import subprocess
import argparse
import tempfile as tmp
import momi
from operator import itemgetter
#import pickle
import dill as pickle
print('Libraries loaded...')


#### MODEL DEFINITIONS ####
def twopop_split(sfs, pop_list, main_pop):
	if len(pop_list) != 2:
		sys.exit('pop_list does not contain 2 populations for 2pop model')
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	new_pop = tp[0]
	TWOSPLIT = momi.DemographicModel(N_e=1e5)
	TWOSPLIT.set_data(sfs)
	TWOSPLIT.add_time_param("tdiv", upper = 1e6)
	TWOSPLIT.add_leaf(main_pop)
	TWOSPLIT.add_leaf(new_pop)
	TWOSPLIT.move_lineages(new_pop, main_pop, t="tdiv")
	TWOSPLIT.optimize(method="TNC")
	return(TWOSPLIT)

def twopop_split_sizes(sfs, pop_list, main_pop):
	if len(pop_list) != 2:
		sys.exit('pop_list does not contain 2 populations for 2pop model')
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	new_pop = tp[0]
	TWOSIZES = momi.DemographicModel(N_e=1e5)
	TWOSIZES.set_data(sfs)
	TWOSIZES.add_size_param("n_" + main_pop, upper = 1e7)
	TWOSIZES.add_size_param("n_" + new_pop, upper = 1e5)
	TWOSIZES.add_time_param("tdiv", upper = 1e6)
	TWOSIZES.add_leaf(main_pop, N = "n_" + main_pop)
	TWOSIZES.add_leaf(new_pop, N = "n_" + new_pop)
	TWOSIZES.move_lineages(new_pop, main_pop, t="tdiv")
	TWOSIZES.optimize(method="TNC")
	return(TWOSIZES)

def twopop_migration(sfs, pop_list, main_pop):
	if len(pop_list) != 2:
		sys.exit('pop_list does not contain 2 populations for 2pop model')
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	new_pop = tp[0]
	TWOMIGRATIONEVS = momi.DemographicModel(N_e=1e5)
	TWOMIGRATIONEVS.set_data(sfs)
	##
	TWOMIGRATIONEVS.add_time_param("tdiv", upper = 1e6)
	TWOMIGRATIONEVS.add_time_param("tmig_" + main_pop + "_" + new_pop, upper_constraints = ['tdiv'])
	##
	TWOMIGRATIONEVS.add_size_param("n_" + main_pop, upper = 1e7)
	TWOMIGRATIONEVS.add_size_param("n_" + new_pop, upper = 1e5)
	##
	TWOMIGRATIONEVS.add_leaf(main_pop, N = "n_" + main_pop)
	TWOMIGRATIONEVS.add_leaf(new_pop, N = "n_" + new_pop)
	##
	TWOMIGRATIONEVS.add_pulse_param("mfrac_" + main_pop + "_" + new_pop, upper = 0.5)
	##
	TWOMIGRATIONEVS.move_lineages(new_pop, main_pop, t = "tdiv")
	TWOMIGRATIONEVS.move_lineages(new_pop, main_pop, t = "tmig_" + main_pop + "_" + new_pop, p = "mfrac_" + main_pop + "_" + new_pop)
	##
	TWOMIGRATIONEVS.optimize(method="TNC")
	return(TWOMIGRATIONEVS)


def twopop_double_migration(sfs, pop_list, main_pop):
	if len(pop_list) != 2:
		sys.exit('pop_list does not contain 2 populations for 2pop model')
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	new_pop = tp[0]
	TWOMIGRATIONEVS = momi.DemographicModel(N_e=1e5)
	TWOMIGRATIONEVS.set_data(sfs)
	##
	TWOMIGRATIONEVS.add_time_param("tdiv", upper = 1e6)
	TWOMIGRATIONEVS.add_time_param("tmig_" + main_pop + "_" + new_pop + "_A", upper_constraints = ['tdiv'])
	TWOMIGRATIONEVS.add_time_param("tmig_" + main_pop + "_" + new_pop + "_B", upper_constraints = ['tdiv'])
	##
	TWOMIGRATIONEVS.add_size_param("n_" + main_pop, upper = 1e7)
	TWOMIGRATIONEVS.add_size_param("n_" + new_pop, upper = 1e5)
	##
	TWOMIGRATIONEVS.add_leaf(main_pop, N = "n_" + main_pop)
	TWOMIGRATIONEVS.add_leaf(new_pop, N = "n_" + new_pop)
	##
	TWOMIGRATIONEVS.add_pulse_param("mfrac_" + main_pop + "_" + new_pop + "_A", upper = 0.5)
	TWOMIGRATIONEVS.add_pulse_param("mfrac_" + main_pop + "_" + new_pop + "_B", upper = 0.5)
	##
	TWOMIGRATIONEVS.move_lineages(new_pop, main_pop, t = "tdiv")
	TWOMIGRATIONEVS.move_lineages(new_pop, main_pop, t = "tmig_" + main_pop + "_" + new_pop + "_A", p = "mfrac_" + main_pop + "_" + new_pop + "_A")
	TWOMIGRATIONEVS.move_lineages(new_pop, main_pop, t = "tmig_" + main_pop + "_" + new_pop + "_B", p = "mfrac_" + main_pop + "_" + new_pop + "_B")
	##
	TWOMIGRATIONEVS.optimize(method="TNC")
	return(TWOMIGRATIONEVS)


def threepop_split_from_main(sfs, pop_list, main_pop):
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	pair1 = main_pop + "_" + tp[0]
	pair2 = main_pop + "_" + tp[1]
	THREEPSFM = momi.DemographicModel(N_e=1e5)
	THREEPSFM.set_data(sfs)
	THREEPSFM.add_time_param("tdiv_" + pair1)
	THREEPSFM.add_time_param("tdiv_" + pair2)
	THREEPSFM.add_leaf(main_pop)
	THREEPSFM.add_leaf(tp[0])
	THREEPSFM.add_leaf(tp[1])
	THREEPSFM.move_lineages(tp[0],main_pop, t="tdiv_" + pair1)
	THREEPSFM.move_lineages(tp[1],main_pop, t="tdiv_" + pair2)
	THREEPSFM.optimize(method="TNC")
	return(THREEPSFM)

def threepop_w_migration_in_beach(sfs,pop_list,main_pop):
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	pair1 = main_pop + "_" + tp[0]
	pair2 = main_pop + "_" + tp[1]
	THREEMIGBEACH = momi.DemographicModel(N_e=1e5)
	THREEMIGBEACH.set_data(sfs)
	THREEMIGBEACH.add_time_param("tdiv_" + pair1)
	THREEMIGBEACH.add_time_param("tdiv_" + pair2)
	THREEMIGBEACH.add_time_param("tmig_" + tp[0] + "_" + tp[1],upper_constraints = ["tdiv_" + pair1,"tdiv_" + pair2],upper=1e5)
	THREEMIGBEACH.add_leaf(main_pop)
	THREEMIGBEACH.add_leaf(tp[0])
	THREEMIGBEACH.add_leaf(tp[1])
	THREEMIGBEACH.add_size_param("n_" + main_pop,1e7)
	THREEMIGBEACH.add_size_param("n_" + tp[0],1e5)
	THREEMIGBEACH.add_size_param("n_" + tp[1],1e5)
	THREEMIGBEACH.add_pulse_param("mfrac_" + tp[0] + "_" + tp[1])
	THREEMIGBEACH.move_lineages(tp[0],main_pop, t="tdiv_" + pair1)
	THREEMIGBEACH.move_lineages(tp[1], main_pop,t="tdiv_" + pair1)
	THREEMIGBEACH.move_lineages(tp[0],tp[1],t="tmig_" + tp[0] + "_" + tp[1],p="mfrac_" + tp[0] + "_" + tp[1])
	THREEMIGBEACH.optimize(method="TNC")
	return(THREEMIGBEACH)


def threepop_w_migration_main(sfs,pop_list,main_pop):
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	pair1 = main_pop + "_" + tp[0]
	pair2 = main_pop + "_" + tp[1]
	THREEMIGMAIN = momi.DemographicModel(N_e=1e5)
	THREEMIGMAIN.set_data(sfs)
	THREEMIGMAIN.add_time_param("tdiv_" + pair1)
	THREEMIGMAIN.add_time_param("tdiv_" + pair2)
	THREEMIGMAIN.add_time_param("tmig_" + main_pop + "_" + tp[0],upper_constraints = ["tdiv_" + pair1,"tdiv_" + pair2],upper=1e5)
	THREEMIGMAIN.add_time_param("tmig_" + main_pop + "_" + tp[1],upper_constraints = ["tdiv_" + pair1,"tdiv_" + pair2],upper=1e5)
	THREEMIGMAIN.add_leaf(main_pop)
	THREEMIGMAIN.add_leaf(tp[0])
	THREEMIGMAIN.add_leaf(tp[1])
	THREEMIGMAIN.add_size_param("n_" + main_pop,1e7)
	THREEMIGMAIN.add_size_param("n_" + tp[0],1e5)
	THREEMIGMAIN.add_size_param("n_" + tp[1],1e5)
	THREEMIGMAIN.add_pulse_param("mfrac_" + main_pop + "_" + tp[0])
	THREEMIGMAIN.add_pulse_param("mfrac_" + main_pop + "_" + tp[1])
	THREEMIGMAIN.move_lineages(tp[0],main_pop, t="tdiv_" + pair1)
	THREEMIGMAIN.move_lineages(tp[1], main_pop,t="tdiv_" + pair1)
	THREEMIGMAIN.move_lineages(tp[0],main_pop,t="tmig_" + main_pop + "_" + tp[0],p="mfrac_" + main_pop + "_" + tp[0])
	THREEMIGMAIN.move_lineages(tp[1],main_pop,t="tmig_" + main_pop + "_" + tp[1],p="mfrac_" + main_pop + "_" + tp[1])
	THREEMIGMAIN.optimize(method="TNC")
	return(THREEMIGMAIN)

def threepop_sizes(sfs, pop_list, main_pop):
	tp = pop_list
	tp = [i for i in tp if i != main_pop]
	pair1 = main_pop + "_" + tp[0]
	pair2 = main_pop + "_" + tp[1]
	THREESIZES = momi.DemographicModel(N_e=1e5)
	THREESIZES.set_data(sfs)
	THREESIZES.add_time_param("tdiv_" + pair1)
	THREESIZES.add_time_param("tdiv_" + pair2)
	THREESIZES.add_leaf(main_pop)
	THREESIZES.add_leaf(tp[0])
	THREESIZES.add_leaf(tp[1])
	THREESIZES.add_size_param("n_" + main_pop,1e7)
	THREESIZES.add_size_param("n_" + tp[0],1e5)
	THREESIZES.add_size_param("n_" + tp[1],1e5)
	THREESIZES.move_lineages(tp[0],main_pop, t="tdiv_" + pair1)
	THREESIZES.move_lineages(tp[1],main_pop, t="tdiv_" + pair2)
	THREESIZES.optimize(method="TNC")
	return(THREESIZES)


#### SCRIPT ####

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--sfs", help="input momi-style sfs")
	parser.add_argument("--populations", help="comma separated list of 2-3 populations to analyze (ex: POALB,LO,PKBM)")
	parser.add_argument("--prefix", help="prefix for output file", default = 'momi')
	parser.add_argument("--resample", help="Whether to resample sfs prior to model run (good for bootstrap replicates)", default=False)
	parser.add_argument("--main_pop", help="focal population from which other ones split, for example POALB")
	parser.add_argument("--runs", help="number of runs for parameter search, default = 50", default = 50)
	parser.add_argument("--models", help="comma separated list of models to use", default = 'all')
	
	print("Parsing arguments...")

	# Parse arguments
	args = parser.parse_args()
	sfs = args.sfs
	pops = args.populations
	prefix = args.prefix
	resample = args.resample
	main_pop = args.main_pop
	runs = int(args.runs)
	models = args.models

	# Run

	## Setup SFS
	print('Loading sfs...')
	sfs = momi.Sfs.load(sfs)
	if resample:
		print('Resampling sfs...')
		sfs = sfs.resample()

	## Setup pops
	pops = pops.split(',')
	print(pops)

	## Setup models
	if models.lower() != 'all':
		print('Run will include models: %s' % models)
		models = models.split(',')
	else:
		if (len(pops) == 2):
			models = ['twopop_split','twopop_split_sizes','twopop_migration','twopop_double_migration']
		elif (len(pops)  == 3):
			models = ['threepop_split_from_main','threepop_w_migration_in_beach','threepop_w_migration_main','threepop_sizes']
		else:
			exit('Population string does not include 2 or 3 populations')

	## Run
	iterations = range(1,runs+1)
	multi_results = []
	for model in models:
		model_results = []
		for i in iterations:
			print('Processing %s: iteration %s' % (model,i))
			try:
				res = globals()[model](sfs, pops, main_pop)
				model_results.append(res)
			except ValueError:
				print('\tparameter search failed for some reason, perhaps due to timing conflict...moving on')
				continue
		model_Lhood = [s.log_likelihood() for s in model_results]
		best_run = model_results[max(enumerate(model_Lhood), key=itemgetter(1))[0]]
		multi_results.append(best_run)

		#TWOMIGRATIONEVS_L = [s.log_likelihood() for s in TWOMIGRATIONEVS]
		#best_TWOMIGRATIONEVS = TWOMIGRATIONEVS[max(enumerate(TWOMIGRATIONEVS_L), key=itemgetter(1))[0]]

	## Write out results
	with open(prefix + '.dat', 'wb') as f:
		pickle.dump(multi_results, f)


if __name__ == "__main__":
	main()