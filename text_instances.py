import os

main_experiment = {
    "method_list": ['mult', 'pdf', 'cdf', 'simulate_1000', 'simulate_100', 'full'],
    "J_list": [100],  # personas
    "M_list": [50],  # mesas
    "G_list": [2, 3, 4],  # grupos
    "I_list": [2, 3, 4, 5, 10],  # candidatos
    "lambda_list": [0.5],
    "seed_list": [i + 1 for i in range(20)],
    "pinit_list": [-1],
    "experiment_name": "main_experiment"
}

simulate_runs = {
    "method_list": ['simulate_1000', 'simulate_100'],
    "J_list": [100],  # personas
    "M_list": [50],  # mesas
    "G_list": [2, 3, 4],  # grupos
    "I_list": [2, 3, 4, 5, 10],  # candidatos
    "lambda_list": [0.5],
    "seed_list": [i + 1 for i in range(20)],
    "pinit_list": [-1],
    "experiment_name": "simulate_runs"
}

convergence_experiment = {
    "method_list": ['full'],    
    "J_list": [100],  # personas
    "M_list": [50],  # mesas
    "G_list": [2, 3],  # grupos
    "I_list": [2, 3],  # candidatos
    "lambda_list": [0.5],
    "seed_list": [1],
    "pinit_list": [i + 1 for i in range(20)],
    "experiment_name": "convergence_experiment"
}

lambda_experiment = {
    "method_list": ['full'],
    "J_list": [100],  # personas
    "M_list": [50],  # mesas
    "G_list": [2, 3],  # grupos
    "I_list": [2, 3],  # candidatos
    "lambda_list": [0.05 * i for i in range(21)],  # as percentage
    "seed_list": [i + 1 for i in range(20)],
    "pinit_list": [-1],
    "experiment_name": "lambda_experiment"
}

text_experiment = {
    "method_list": ['full'],
    "J_list": [100],  # personas
    "M_list": [50],  # mesas
    "G_list": [2, 3],  # grupos
    "I_list": [2],  # candidatos
    "lambda_list": [0.5],
    "seed_list": [1],
    "pinit_list": [1],
    "experiment_name": "text_experiment"
}

# Function to generate run instances
def print_run_instances(experiment_params, file_name=None):
    if file_name is None:
        file_name = experiment_params['experiment_name'] + '.txt'

    # Ensure the script directory exists
    os.makedirs('experiment_scripts', exist_ok=True)
    
    file_path = os.path.join('experiment_scripts', file_name)

    with open(file_path, 'w') as f:
        for method in experiment_params["method_list"]:
            f.write(f'\n')
            for J in experiment_params["J_list"]:
                for M in experiment_params["M_list"]:
                    for G in experiment_params["G_list"]:
                        for I in experiment_params["I_list"]:
                            for lambda_ in experiment_params["lambda_list"]:
                                for seed in experiment_params["seed_list"]:
                                    for pinit in experiment_params["pinit_list"]:
                                        f.write(f'python3 EM_cluster.py {method} {J:d} {M:d} {G:d} {I:d} {lambda_:.2f} {seed:d} {pinit:d}\n')
                                    
    print(f"Run instances for {experiment_params['experiment_name']} saved to {file_path}")


# Example usage
if __name__ == "__main__":
    # Select which experiment to run and generate the corresponding script
    print("Generating scripts for different experiments...")

    # Generate main experiment script
    print_run_instances(main_experiment)

    # Generate convergence experiment script
    print_run_instances(convergence_experiment)

    # Generate lambda experiment script
    print_run_instances(lambda_experiment)

    # Generate simulate runs script
    print_run_instances(simulate_runs)

    print("All scripts generated and saved in 'experiment_scripts' directory.")
