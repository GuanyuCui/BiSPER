# BiSPER

*(This is the official GitHub repository for the KDD 2025 paper: Mixing Time Matters: Accelerating Effective Resistance Estimation via Bidirectional Method.)*


Welcome to the official BiSPER documentation. This guide provides a comprehensive walkthrough to replicate our experiments, covering everything from dataset preparation and code compilation to running experiments and plotting results.

## Environment Setup

### Requirements
- **Tools**: gcc (version 7.5 or newer), Python 3.x, Eigen 3.4

### Setup Instructions
1. Clone the repository to your local machine.
2. Install the necessary Python packages using the command below:

    ```shell
    pip install numpy scipy pandas matplotlib seaborn networkx
    ```

## Dataset Preparation

### Preparing Datasets

**Note**: To save time, it is recommended to proceed directly to step 3.

1. **Original Real-World Datasets**: Obtain these from the [SNAP dataset page](https://snap.stanford.edu/data/). The initial run of our code will preprocess these datasets, generating files named `<dataset_name>-compressed_sorted.bin`.

2. **Generated Synthetic Datasets**: Create Erdos-Renyi random graphs by executing `generate_synthetic.py`. This process generates an edge list file named `synthetic.txt`.

3. **Preprocessed Datasets**: For convenience, you can directly download our preprocessed datasets (in adjacency list and binary formats) from [this link](https://mega.nz/folder/pPghxKpZ#ZZCDT2H844otXKrchH2jbA). 

Ensure all datasets are placed in the `BiSPER/datasets` directory.

**Credit**: Original datasets are courtesy of [SNAP](https://snap.stanford.edu/data/).

### Calculating Lambda Values (Optional)

Pre-calculated lambda values are found in `BiSPER/datasets/<dataset_name>.lambda`. To recalculate:

- Ensure the existence of `<dataset_name>-compressed_sorted.bin` (refer to Dataset Preparation).
- Execute `python get_lambda.py <dataset_name>` from the BiSPER root directory. The lambda value will be recalculated and saved to `BiSPER/datasets/<dataset_name>.lambda`.

## Compilation

Navigate to the `BiSPER/` directory and compile the code with the following commands:

```shell
cd src
g++ -std=c++14 -O3 *.cpp -o SPER -lstdc++fs -pthread
mv SPER ..
cd ..
```

**Credit**: Code for Eigen 3.4 is taken from [the official website](https://eigen.tuxfamily.org/) of the Eigen library, code for AMC / GEER algorithms are from [their official GitHub page](https://github.com/AnryYang/GEER), code for Bipush / Push / AbWalk algorithms are from [their official GitHub page](https://github.com/mhliao516/Resistance-Landmark), and code for Bipush-vl / Push-vl / RW-vl algorithms are from [their official GitHub page](https://github.com/mhliao0516/EffectiveResistanceMultipleLandmark).

## Running Experiments

### Reproducing Our Experiments

Follow these steps to replicate our experiments:

1. Grant execution permissions to the shell scripts:

    ```shell 
    chmod +x *.sh
    ```

2. Execute a specific experiment script:

    **Note**: Experiment 0 stands for the $L_{\max}$ truncated value v.s. $L_{\max}$ experiment in Table 2.

    ```shell
    ./Experiment-[0/I/II/III].sh
    ```

### Custom Experiments with SPER

You can utilize the compiled SPER executable for custom experiments. The usage syntax is:

```shell
SPER 
    [--dataset <dataset_name>] 
    [--algorithm <algorithm_name>] 
    [--num_query <num_query>] 
    [--L_max <L_max>] 
    [--eps <eps>] 
    [--p_f <p_f>]
    [--num_landmarks <num_landmarks>]
    [--num_samples <num_samples>] 
    [--r_max <r_max>]
```

#### Argument Details

Below are detailed explanations of the arguments that can be used with the program:

- `--dataset`: Specifies the dataset to be used. The program initially attempts to locate the compressed and sorted binary file at `BiSPER/datasets/<dataset_name>-compressed_sorted.bin`. If this file cannot be found, it then searches for a text file at `BiSPER/datasets/<dataset_name>.txt`. Should both files be unavailable, the program will raise an exception. The default dataset is `Facebook`.

- `--algorithm`: Determines the algorithm to be employed. Currently supported algorithms include: `BiSPER`, `AMC`, `GEER`, `Bipush`, and `Push`. The default algorithm is `BiSPER`.

- `--num_query`: Sets the number of query pairs to generate. Both the pre-sampled queries and their corresponding ground-truth values are stored in the `BiSPER/samples/` directory. The default number is `100`.

- `--L_max`: Defines the maximum number of steps, denoted as $L_{\max}$. Users can input an integer or specify `auto`, which calculates $L_{\max}$ based on lambda and epsilon values. This argument is applicable only to the `BiSPER`, `AMC`, and `GEER` algorithms. The default value is `100`.

- `--eps`: Sets the desired absolute error guarantee, $\epsilon$. This parameter is relevant only for the `BiSPER`, `AMC`, and `GEER` algorithms. The default value is ```1e-2```.

- `--num_landmarks`: Indicates the number of landmark nodes in the `Bipush-vl`, `Push-vl`, and `RW-vl` algorithms. The default value is `100`.

- `--num_samples`: Indicates the number of random walks to sample in the `Bipush` algorithm. The default value is `10,000`.

- `--r_max`: Specifies the push threshold for both the `Bipush` and `Push` algorithms. The default value is ```1e-4```.

Results are stored in `BiSPER/results/` as `.out` files for any text editor access.

## Plotting Line Charts (Optional)

To recreate the line charts in experiments from our papers, execute the following command within the `BiSPER/results/` directory:

```shell
python plot.py
```

This generates and saves PGF files for visualizations.