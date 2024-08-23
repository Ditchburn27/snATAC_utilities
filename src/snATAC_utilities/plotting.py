#### Plotting Functions for snATAC-seq 
#####################################
## Plot cell counts of specified obs for 
## different tsse thresholds. 
def plot_tsse_countplots(adata, groupby_observation, batch_key=None, log_scale=False):
    """
    Plots cell counts of specified obs for different tsse thresholds,
    optionally further grouped by 'batch' or a custom key.
    
    Parameters:
    - adata: Anndata object containing the observations data.
    - groupby_observation: The key in adata.obs to group by (e.g., 'cell_type').
    - batch_key: Optional; The key in adata.obs to use for further grouping (default: None).
    - log_scale: Boolean, if True, sets y-axis to log scale.
    """
    if batch_key is None:
        _plot_tsse_countplots_single(adata, groupby_observation, log_scale)
    else:
        # Unique batches in the data
        batches = adata.obs[batch_key].unique()
        num_batches = len(batches)
        num_plots = 4  # We always want 4 plots per batch
        
        # Create figure and subplots based on number of batches and num_plots per row
        fig, axes = plt.subplots(num_batches, num_plots, figsize=(16, 3 * num_batches), sharey=True)
        
        # Define the conditions for subplots
        conditions = [
            (adata.obs['tsse'] < 5),
            (adata.obs['tsse'] >= 5) & (adata.obs['tsse'] < 10),
            (adata.obs['tsse'] >= 10) & (adata.obs['tsse'] < 20),
            (adata.obs['tsse'] >= 20)
        ]
        
        # Titles for subplots
        titles = [
            'tsse < 5',
            '5 <= tsse < 10',
            '10 <= tsse < 20',
            'tsse >= 20'
        ]
        
        # Loop over batches and plot for each batch
        for i, batch in enumerate(batches):
            # Loop over conditions (tsse thresholds)
            for j in range(num_plots):
                ax = axes[i, j]  # Select the correct subplot
                
                cond = conditions[j]
                
                sns.countplot(x=groupby_observation, data=adata.obs[cond & (adata.obs[batch_key] == batch)],
                              ax=ax, palette='Set3')
                
                ax.set_title(f'TSSE {titles[j]}, {batch_key}: {batch}')
                ax.set_xlabel(groupby_observation)
                ax.set_ylabel('Count')
                if log_scale:
                    ax.set_yscale('log')
                plt.setp(ax.get_xticklabels(), rotation=45, ha='right')  # Rotate x-axis labels
        
        # Adjust layout and show plot
        plt.tight_layout()
        plt.show()

def _plot_tsse_countplots_single(adata, groupby_observation, log_scale=False):
    """
    Helper function to plot cell counts of specified obs for different tsse thresholds
    without further grouping by batch.
    
    Parameters:
    - adata: Anndata object containing the observations data.
    - groupby_observation: The key in adata.obs to group by (e.g., 'cell_type').
    - log_scale: Boolean, if True, sets y-axis to log scale.
    """
    # Create figure and 2x2 subplots with smaller figsize
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharey=True)
    axes = axes.flatten()
    
    # Define the conditions and titles for subplots
    conditions = [
        (adata.obs['tsse'] < 5),
        (adata.obs['tsse'] >= 5) & (adata.obs['tsse'] < 10),
        (adata.obs['tsse'] >= 10) & (adata.obs['tsse'] < 20),
        (adata.obs['tsse'] >= 20)
    ]
    titles = [
        'tsse < 5',
        '5 <= tsse < 10',
        '10 <= tsse < 20',
        'tsse >= 20'
    ]
    
    for i, ax in enumerate(axes):
        cond = conditions[i]
        sns.countplot(x=groupby_observation, data=adata.obs[cond], ax=ax, palette='Set3')  # Set different colors
        ax.set_title(titles[i])
        ax.set_xlabel(groupby_observation)
        ax.set_ylabel('Count')
        if log_scale:
            ax.set_yscale('log')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')  # Rotate x-axis labels
    
    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()

#####################################
## Function to plot a grid of subplots for counting cell numbers
## for each label in the provided count key
def plot_grouped_counts(adata, group_key, count_key, num_cols=5):
    """
    Plots subplots of count plots for each unique group in adata.obs.

    Parameters:
    - adata: Anndata object containing the observations data.
    - group_key: The key in adata.obs to group by (e.g., 'batch').
    - count_key: The key in adata.obs to count (e.g., 'major_cluster').
    - num_cols: Number of columns in the subplot grid.
    """
    unique_groups = adata.obs[group_key].unique()
    num_subplots = len(unique_groups)
    num_rows = -(-num_subplots // num_cols)  # Ceiling division to ensure at least num_subplots number of rows
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5*num_cols, 5*num_rows), sharey=True)
    # Flatten the axes array if it's not already 1-dimensional
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    # Define a color palette for each unique count_key
    palette = sns.color_palette("husl", len(adata.obs[count_key].unique()))
    # Iterate over each unique group
    for i, group in enumerate(unique_groups):
        # Filter data for current group
        group_data = adata.obs[adata.obs[group_key] == group]    
        # Count plot for each count_key label
        sns.countplot(x=count_key, data=group_data, ax=axes[i], palette=palette)
        axes[i].set_title(f'{group_key.capitalize()} {group}')      
        # Rotate x-labels for better readability
        axes[i].tick_params(axis='x', rotation=45)        
        # Set labels for better readability
        axes[i].set_xlabel(count_key.capitalize())
        axes[i].set_ylabel('Count')
    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.show()

#####################################
# Function for plotting scatter plot
# of adata.obs columns
def scatter_plot_by_group(adata, x_var, y_var, groupby_var=None):
    """
    Plot scatter plots of x_var vs y_var grouped by groupby_var if provided.

    Parameters:
    adata (AnnData): Annotated data object containing the data.
    x_var (str): Name of the observation column for the x-axis.
    y_var (str): Name of the observation column for the y-axis.
    groupby_var (str or None): Name of the observation column for grouping.
                               If None, plot a single plot for the whole dataset.
    """
    if groupby_var is None:
        # Plot a single plot for the whole dataset
        fig, ax = plt.subplots(figsize=(8, 4))
        x = adata.obs[x_var]
        y = adata.obs[y_var]

        # Calculate the trend line
        coefficients = np.polyfit(x, y, 1)
        polynomial = np.poly1d(coefficients)
        trend_line = polynomial(x)

        # Calculate the correlation coefficient
        correlation_coefficient, p_value = pearsonr(x, y)

        # Create scatter plot
        ax.scatter(x, y, alpha=0.2)
        ax.plot(x, trend_line, color='grey', label=f'Correlation (r={correlation_coefficient:.2f})')

        # Add titles and labels
        ax.set_title(f'{x_var} vs {y_var}', fontsize=14)
        ax.set_xlabel(x_var, fontsize=12)
        ax.set_ylabel(y_var, fontsize=12)

        # Customize tick labels
        ax.tick_params(axis='both', labelsize=10)

        # Add legend
        ax.legend(fontsize=10)

        # Show plot
        plt.show()

    else:
        # Group by groupby_var and plot separate subplots
        categories = adata.obs[groupby_var].cat.categories
        n_rows = len(categories) // 2 + len(categories) % 2
        n_cols = 2

        fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        axs = axs.flatten()

        # Iterate over each category
        for i, category in enumerate(categories):
            ax = axs[i]
            subset = adata[adata.obs[groupby_var] == category]
            x = subset.obs[x_var]
            y = subset.obs[y_var]

            # Calculate the trend line
            coefficients = np.polyfit(x, y, 1)
            polynomial = np.poly1d(coefficients)
            trend_line = polynomial(x)

            # Calculate the correlation coefficient
            correlation_coefficient, p_value = pearsonr(x, y)

            # Get color for the category from adata.uns[groupby_var + '_colors'] if available
            if f'{groupby_var}_colors' in adata.uns:
                colors = adata.uns[f'{groupby_var}_colors']
                color = colors[i % len(colors)]
            else:
                color = 'blue'  # Default color if colors not available

            # Create scatter plot with the corresponding color
            ax.scatter(x, y, alpha=0.2, color=color)
            ax.plot(x, trend_line, color='grey', label=f'Correlation (r={correlation_coefficient:.2f})')

            # Add titles and labels
            ax.set_title(f'{x_var} vs {y_var} - {groupby_var}: {category}', fontsize=16)
            ax.set_xlabel(x_var, fontsize=12)
            ax.set_ylabel(y_var, fontsize=12)

            # Customize tick labels
            ax.tick_params(axis='both', labelsize=10)

            # Add legend
            ax.legend(fontsize=10)

        # Remove any extra empty subplots
        for i in range(len(categories), len(axs)):
            fig.delaxes(axs[i])

        # Adjust layout
        plt.tight_layout()
        plt.subplots_adjust(top=0.93)

        # Add main figure title
        plt.suptitle(f'Scatter plot of {x_var} vs {y_var} by {groupby_var}', fontsize=20)

        # Show plot
        plt.show()

#####################################
# Function for plotting scatter plot
# of adata.obs columns with a list 
# of columns for the y variable
def scatter_plot_multiple_y(adata, x_var, y_vars):
    """
    Plot scatter plots of x_var vs each y_var in y_vars with data points colored differently for each subplot.

    Parameters:
    adata (AnnData): Annotated data object containing the data.
    x_var (str): Name of the observation column for the x-axis.
    y_vars (list of str): List of observation column names for the y-axis.
    """
    n_plots = len(y_vars)
    n_rows = n_plots // 2 + n_plots % 2
    n_cols = 2

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    axs = axs.flatten()

    # Get a colormap for the plots
    cmap = plt.cm.get_cmap('Set2')

    # Loop through each y variable and assign different colors for each subplot
    for i, y_var in enumerate(y_vars):
        ax = axs[i]
        x = adata.obs[x_var]
        y = adata.obs[y_var]

        # Calculate the trend line
        coefficients = np.polyfit(x, y, 1)
        polynomial = np.poly1d(coefficients)
        trend_line = polynomial(x)

        # Calculate the correlation coefficient
        correlation_coefficient, p_value = pearsonr(x, y)

        # Plot scatter plot with different colors for each subplot and add transparency
        color = cmap(i / n_plots)  # Assign color based on subplot index
        ax.scatter(x, y, color=color, alpha=0.2, edgecolor='none')

        # Plot trend line
        ax.plot(x, trend_line, color='grey', label=f'Correlation (r={correlation_coefficient:.2f})')

        # Add titles and labels
        ax.set_title(f'{x_var} vs {y_var}', fontsize=14)
        ax.set_xlabel(x_var, fontsize=12)
        ax.set_ylabel(y_var, fontsize=12)

        # Customize tick labels
        ax.tick_params(axis='both', labelsize=10)

        # Add legend
        ax.legend(fontsize=10)

    # Remove any extra empty subplots
    for i in range(n_plots, len(axs)):
        fig.delaxes(axs[i])

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    # Add main figure title outside subplots
    plt.suptitle(f'Scatter plot of {x_var} vs multiple observations', fontsize=20)

    # Show plot
    plt.show()