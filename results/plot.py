import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

def data2dict(data_str):
    pairs = data_str.split()
    parsed_dict = {}
    for pair in pairs:
        key, value = pair.split(':')
        parsed_dict[key] = value
        
    return parsed_dict

def plot_experiment_I(x, y, L_max : str, algorithm2marker, algorithm2color):
    dataset_names = ['Facebook', 'DBLP', 'Youtube', 'Orkut', 'LiveJournal', 'Friendster']
    algorithm_names = ['BiSPER', 'GEER', 'AMC']

    columns = ['Algorithm', 'L_max', 'eps', 'max.error', 'avg.error', 'avg.time(ms)']
    numeric_columns = ['eps', 'max.error', 'avg.error', 'avg.time(ms)']
    
    sns.set_theme(style = "ticks")
    num_subplot_cols = 3
    num_subplot_rows = 2 # math.ceil(len(dataset_names) / num_subplot_cols)
    fig, axes = plt.subplots(num_subplot_rows, num_subplot_cols, figsize = (16, 8))
    
    # List to store added algorithm labels.
    added_labels = []
    
    if y not in ['max.error', 'avg.error', 'eps']:
        raise
    if y == 'max.error':
        y_label_name = 'Maximum Absolute Error'
    elif y == 'avg.error':
        y_label_name = 'Average Absolute Error'
    elif y == 'eps':
        y_label_name = r'Epsilon ($\epsilon$)'
        
    if x not in ['avg.time(ms)', 'max.error', 'avg.error']:
        raise
    if x == 'avg.time(ms)':
        x_label_name = 'Average Query Time (ms)'
    elif x == 'max.error':
        x_label_name = 'Maximum Absolute Error'
    elif x == 'avg.error':
        x_label_name = 'Average Absolute Error'
    
    if x == 'avg.time(ms)':
        x_lims = [(1e1, 1e4), (1e1, 1e6), (1e1, 1e6), (1e2, 1e5), (1e1, 1e6), (1e2, 1e6)]
    if y == 'max.error':
        y_lims = [(1e-5, 1e-1), (1e-5, 1e-1), (1e-5, 1e-0), (1e-7, 1e-1), (1e-5, 1e-0), (1e-6, 1e-0)]
    elif y == 'avg.error':
        y_lims = [(1e-5, 1e-2), (1e-6, 1e-1), (1e-6, 1e-1), (1e-7, 1e-2), (1e-6, 1e-1), (1e-7, 1e-1)]
    
    for (i, dataset_name) in zip(range(len(dataset_names)), dataset_names):
        # Open output file.
        f = open(f'{dataset_name}.out')
        # Read contents.
        data = f.readlines()
        # To list, to DataFrame.
        dict_list = [data2dict(data_str) for data_str in data]
        # Filter columns.
        df = pd.DataFrame(dict_list)[columns]
        df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors = 'coerce')
        # Filter data in Experiment I.
        df = df[df['L_max'] == L_max]
                
        for algorithm_name in algorithm_names:
            tmp_df = df[df['Algorithm'] == algorithm_name]
            tmp_df = tmp_df.sort_values(by = 'eps')
            ax = axes[i // num_subplot_cols, i % num_subplot_cols]
            # Check if label is already added.
            if algorithm_name not in added_labels: 
                ax.plot(tmp_df[x], tmp_df[y], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', label = algorithm_name, color = algorithm2color[algorithm_name])
                # Add label to added_labels list.
                added_labels.append(algorithm_name)
            else:
                ax.plot(tmp_df[x], tmp_df[y], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', color = algorithm2color[algorithm_name])
        
        if y == 'eps' and (x == 'max.error' or x == 'avg.error'):
            ax.plot([10 ** (-i) for i in range(1, 4)], [10 ** (-i) for i in range(1, 4)] , linestyle = '--', color = 'gray')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.tick_params(axis = 'x', labelsize = 14)
        ax.tick_params(axis = 'y', labelsize = 14)
        
        ax.set_xlabel(x_label_name, fontsize = 14)

        # Set y-axis label only for the first column.
        if i % num_subplot_cols == 0:
            ax.set_ylabel(y_label_name, fontsize = 14)
    
        # Set x/y-axis limits.
        if L_max == '100':
            if y == 'max.error' or y == 'avg.error':
                ax.set_ylim(y_lims[i][0], y_lims[i][1])
            if x == 'avg.time(ms)':
                ax.set_xlim(x_lims[i][0], x_lims[i][1])
            
        ax.set_title(dataset_name, fontsize = 14)

            
    fig.legend(bbox_to_anchor = (0.5, 0.97), loc = 'upper center', ncols = len(algorithm_names), frameon = False, fontsize = 14)
    # Adjusting vertical and horizontal spacing.
    plt.subplots_adjust(hspace = 0.35, wspace = 0.2)
    plt.savefig('Experiment-I-results-{}-{}.pgf'.format(x, y), bbox_inches = 'tight')
    # plt.show()
    
def plot_experiment_II(algorithm2marker, algorithm2color):
    algorithm_names = ['BiSPER', 'GEER', 'Bipush', 'Push', 'AbWalk', 'Bipush-vl', 'Push-vl', 'RW-vl']
    columns = ['Algorithm', 'L_max', 'num_samples', 'r_max', 'max.error', 'avg.error', 'avg.time(ms)']
    numeric_columns = ['num_samples', 'r_max', 'max.error', 'avg.error', 'avg.time(ms)']
    
    sns.set_theme(style = "ticks")
    fig = plt.figure(figsize = (4, 4))
    ax = fig.add_subplot(111)
    
    # Open output file.
    f = open('Facebook.out')
    # Read contents.
    data = f.readlines()
    # To list, to DataFrame.
    dict_list = [data2dict(data_str) for data_str in data]
    # Filter columns.
    df = pd.DataFrame(dict_list)[columns]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors = 'coerce')
    # Filter data in Experiment II.
    df = df[df['L_max'] == 'auto']
    
    for algorithm_name in algorithm_names:
        if algorithm_name in ['Bipush']:
            r_maxs = ['1e-4', '1e-5', '1e-6']
            for r_max in r_maxs:
                # Get data of some algorithm, and sort values.
                tmp_df = df[(df['Algorithm'] == algorithm_name) & (df['r_max'] == float(r_max))]
                tmp_df = tmp_df.sort_values(by = 'num_samples')
                if r_max == r_maxs[0]:
                    # Plot.
                    ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', label = algorithm_name, color = algorithm2color[algorithm_name])
                else:
                    ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', color = algorithm2color[algorithm_name])

                ax.annotate(fr'$r_{{\max}} = 10^{{{r_max[-2:]}}}$',
                        xy = (tmp_df['avg.time(ms)'].iloc[-1], tmp_df['avg.error'].iloc[-1]),
                        xytext = (12, 12),
                        textcoords = 'offset points',
                        fontsize = 13,
                        usetex = True,
                        arrowprops = dict(arrowstyle = '->', color = 'black'))
                    
        else:
            # Get data of some algorithm, and sort values.
            tmp_df = df[df['Algorithm'] == algorithm_name]
            tmp_df = tmp_df.sort_values(by = 'avg.error')
            # Plot.
            ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', label = algorithm_name, color = algorithm2color[algorithm_name])
        
    ax.set_ylim(1e-11, 1e-2)
    ax.set_xlim(1e3, 1e6)
        
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis = 'x', labelsize = 14)
    ax.tick_params(axis = 'y', labelsize = 14)
    ax.set_ylabel('Average Absolute Error', fontsize = 14)
    ax.set_xlabel('Average Query Time (ms)', fontsize = 14)

    # Get legends' handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    # Row-first
    nrow = 2
    ncol = math.ceil(len(handles) / nrow)

    handles_labels = list(zip(handles, labels))
    handles_labels_sorted = [
        handles_labels[i * ncol + j]
        for j in range(ncol)
        for i in range(nrow)
        if i * ncol + j < len(handles_labels)  # Index out of range
    ]
    handles_sorted, labels_sorted = zip(*handles_labels_sorted)
            
    fig.legend(handles_sorted, labels_sorted, bbox_to_anchor = (0.5, 1.05), loc = 'upper center', ncols = ncol, frameon = False, fontsize = 14)
    # Make it square.
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable = 'box')
    plt.savefig('Experiment-II-results.pgf', bbox_inches = 'tight')
    
def plot_experiment_III(algorithm2marker, algorithm2color):
    algorithm_names = ['BiSPER', 'GEER', 'Bipush', 'Push', 'AbWalk', 'Bipush-vl', 'Push-vl', 'RW-vl']
    columns = ['Algorithm', 'L_max', 'num_samples', 'r_max', 'max.error', 'avg.error', 'avg.time(ms)']
    numeric_columns = ['num_samples', 'r_max', 'max.error', 'avg.error', 'avg.time(ms)']
    
    sns.set_theme(style = "ticks")
    fig = plt.figure(figsize = (4, 4))
    ax = fig.add_subplot(111)
    
    # Open output file.
    f = open('synthetic.out')
    # Read contents.
    data = f.readlines()
    # To list, to DataFrame.
    dict_list = [data2dict(data_str) for data_str in data]
    # Filter columns.
    df = pd.DataFrame(dict_list)[columns]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors = 'coerce')
    # Filter data in Experiment III.
    df = df[df['L_max'] == 'auto']
    
    for algorithm_name in algorithm_names:
        if algorithm_name in ['Bipush']:
            r_maxs = ['1e-4', '1e-5', '1e-6']
            for r_max in r_maxs:
                # Get data of some algorithm, and sort values.
                tmp_df = df[(df['Algorithm'] == algorithm_name) & (df['r_max'] == float(r_max))]
                tmp_df = tmp_df.sort_values(by = 'num_samples')
                if r_max == r_maxs[0]:
                    # Plot.
                    ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', label = algorithm_name, color = algorithm2color[algorithm_name])
                else:
                    ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', color = algorithm2color[algorithm_name])
                    
                ax.annotate(fr'$r_{{\max}} = 10^{{{r_max[-2:]}}}$',
                            xy = (tmp_df['avg.time(ms)'].iloc[-1], tmp_df['avg.error'].iloc[-1]),
                            xytext = (12, 12),
                            textcoords = 'offset points',
                            fontsize = 13,
                            usetex = True,
                            arrowprops = dict(arrowstyle = '->', color = 'black')
                )
        else:
            # Get data of some algorithm, and sort values.
            tmp_df = df[df['Algorithm'] == algorithm_name]
            tmp_df = tmp_df.sort_values(by = 'avg.error')
            # Plot.
            ax.plot(tmp_df['avg.time(ms)'], tmp_df['avg.error'], marker = algorithm2marker[algorithm_name], markersize = 7, markeredgewidth = 1.3, markerfacecolor = 'none', label = algorithm_name, color = algorithm2color[algorithm_name])
            
    ax.set_ylim(1e-11, 1e-2)
    ax.set_xlim(1e-3, 1e5)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis = 'x', labelsize = 14)
    ax.tick_params(axis = 'y', labelsize = 14)
    ax.set_ylabel('Average Absolute Error', fontsize = 14)
    ax.set_xlabel('Average Query Time (ms)', fontsize = 14)

    # Get legends' handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    # Row-first
    nrow = 2
    ncol = math.ceil(len(handles) / nrow)

    handles_labels = list(zip(handles, labels))
    handles_labels_sorted = [
        handles_labels[i * ncol + j]
        for j in range(ncol)
        for i in range(nrow)
        if i * ncol + j < len(handles_labels)  # Index out of range
    ]
    handles_sorted, labels_sorted = zip(*handles_labels_sorted)
            
    fig.legend(handles_sorted, labels_sorted, bbox_to_anchor = (0.5, 1.05), loc = 'upper center', ncols = ncol, frameon = False, fontsize = 14)
    # Make it square.
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable = 'box')
    plt.savefig('Experiment-III-results.pgf', bbox_inches = 'tight')
    
if __name__ == '__main__':    
    algorithm2marker = {
    'BiSPER': 'D',
    'GEER': 's',
    'AMC': 'o',
    'Bipush': 'p',
    'Push': 'h',
    'AbWalk': 'v',
    'Bipush-vl': 'x',
    'Push-vl': '*',
    'RW-vl': '^'
    }

    algorithm2color = {
        'BiSPER': 'red',
        'GEER': 'blue',
        'AMC': 'orange',
        'Bipush': 'purple',
        'Push': 'teal',
        'AbWalk': 'brown',
        'Bipush-vl': 'cyan',
        'Push-vl': 'green',
        'RW-vl': 'magenta'
    }
    plot_experiment_I(x = 'avg.time(ms)', y = 'avg.error', L_max = '100', algorithm2marker = algorithm2marker, algorithm2color = algorithm2color)
    plot_experiment_II(algorithm2marker = algorithm2marker, algorithm2color = algorithm2color)
    plot_experiment_III(algorithm2marker = algorithm2marker, algorithm2color = algorithm2color)