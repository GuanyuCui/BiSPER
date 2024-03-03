import itertools
import subprocess
from argparse import ArgumentParser

if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('--experiment', type = str, default = 'I', help = 'Experiment number. (default = 1)')
	parser.add_argument('--dataset', type = str, default = 'Facebook', help = 'Dataset. (default = Facebook)')
	parser.add_argument('--algorithm', type = str, default = 'Power-Iteration', help = 'Algorithm name. (default = Power-Iteration)')
	parser.add_argument('--num_query', type = str, default = '100', help = 'Number of query pairs. (default = 100)')
	parser.add_argument('--L_max', type = str, default = '100', help = 'L_max. (default = 100)')
	args = parser.parse_args()

	experiment = args.experiment
	dataset = args.dataset
	algorithm = args.algorithm
	num_query = args.num_query
	L_max = args.L_max

	program_path = './SPER' 
	if experiment == 'I':
		assert dataset == 'Facebook' and algorithm == 'Power-Iteration' and num_query == '1'
		command = [program_path, '--dataset', 'Facebook', '--algorithm', 'Power-Iteration', '--num_query', '1', '--L_max', L_max, '--eps', '0.0']
		try:
			subprocess.run(command, check = True)
		except subprocess.CalledProcessError as e:
			print(f"Error calling C++ program: {e}")
	elif experiment == 'II':
		if algorithm == 'AMC':
			epsilons = ['1e-1', '2e-1', '5e-1']
		elif algorithm == 'GEER' and dataset == 'Friendster':
			epsilons = ['1e-2', '2e-2', '5e-2', '1e-1', '2e-1', '5e-1']
		else:
			epsilons = ['1e-3', '2e-3', '5e-3', '1e-2', '2e-2', '5e-2', '1e-1']
		for eps in epsilons:
			command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--L_max', L_max, '--eps', eps]
			try:
				subprocess.run(command, check = True)
			except subprocess.CalledProcessError as e:
				print(f"Error calling C++ program: {e}")
	elif experiment == 'III':
		assert dataset == 'Facebook'
		dataset = 'Facebook'
		if algorithm == 'Bipush':
			num_sampless = ['1000', '10000', '100000']
			r_maxs = ['1e-6', '1e-5', '1e-4']
			for r_max, num_samples in itertools.product(r_maxs, num_sampless):
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--num_samples', num_samples, '--r_max', r_max]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")
		elif algorithm == 'Push':
			r_maxs = ['1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '1e-5', '1e-4']
			for r_max in r_maxs:
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--num_samples', '0', '--r_max', r_max]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")
		else:
			if algorithm == "GEER":
				epsilons = ['1e-3', '1e-2', '1e-1']
			else:
				epsilons = ['1e-7', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1']
			for eps in epsilons:
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--L_max', 'auto', '--eps', eps]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")
	elif experiment == 'IV':
		assert dataset == 'ER'
		dataset = 'ER'
		if algorithm == 'Bipush':
			num_sampless = ['1000', '10000']
			r_maxs = ['1e-6', '1e-5', '1e-4']
			for r_max, num_samples in itertools.product(r_maxs, num_sampless):
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--num_samples', num_samples, '--r_max', r_max]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")
		elif algorithm == 'Push':
			r_maxs = ['1e-7', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1']
			for r_max in r_maxs:
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--num_samples', '0', '--r_max', r_max]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")
		else:
			epsilons = ['1e-7', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1']
			for eps in epsilons:
				command = [program_path, '--dataset', dataset, '--algorithm', algorithm, '--num_query', num_query, '--L_max', 'auto', '--eps', eps]
				try:
					subprocess.run(command, check = True)
				except subprocess.CalledProcessError as e:
					print(f"Error calling C++ program: {e}")