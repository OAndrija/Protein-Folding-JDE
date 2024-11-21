import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('data/results.csv')

# Sort data by the 'run' column
data = data.sort_values(by='run')

#Save the sorted data to a new CSV file
data.to_csv('data/sorted_results.csv', index=False)

# Extract sorted data columns
seeds = data['seed']
E_values = data['E']
runtime_values = data['runtime']
speed_values = data['speed']

# Set y-axis limits
E_ylim = (-5.6, -3)
runtime_ylim = (3, 6.5)
speed_ylim = (150000, 350000)

# Create subplots
fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Function to plot data with stats
def plot_with_stats(ax, x, y, y_label, color, ylim):
    y_min = y.min()
    y_max = y.max()
    y_mean = y.mean()
    y_std = y.std()

    ax.plot(x, y, marker='o', color=color)
    ax.set_ylabel(y_label)
    ax.set_ylim(ylim)

    ax.axhline(y_mean, color=color, linestyle='--', label=f'Average: {y_mean:.2f}')
    ax.fill_between(x, y_mean - y_std, y_mean + y_std, color=color, alpha=0.2, label=f'Std Dev: {y_std:.2f}')

    ax.plot(x[y.idxmin()], y_min, marker='o', color='orange', label=f'Min: {y_min:.2f}')
    ax.plot(x[y.idxmax()], y_max, marker='o', color='yellow', label=f'Max: {y_max:.2f}')
    ax.legend()

# Plot each metric
plot_with_stats(axes[0], seeds, E_values, "E", 'b', E_ylim)
axes[0].set_title("E Values for Each Seed Run")

plot_with_stats(axes[1], seeds, runtime_values, "Runtime (s)", 'g', runtime_ylim)
axes[1].set_title("Runtime for Each Seed Run")

plot_with_stats(axes[2], seeds, speed_values, "Speed", 'r', speed_ylim)
axes[2].set_title("Speed for Each Seed Run")
axes[2].set_xlabel("Seed")

# Adjust x-axis ticks for better readability
axes[2].set_xticks(seeds[::5])

# Final layout adjustments
plt.tight_layout()

# Save the figure
plt.savefig('data/graph_vaja1.png')

# Show the figure
plt.show()
