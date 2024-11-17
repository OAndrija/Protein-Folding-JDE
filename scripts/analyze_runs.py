import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('data/results.csv')

seeds = data['seed']
E_values = data['E']
runtime_values = data['runtime']
speed_values = data['speed ']

E_ylim = (-5.6, -3)
runtime_ylim = (1.8, 2.3)
speed_ylim = (450000, 550000)

fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

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

plot_with_stats(axes[0], seeds, E_values, "E", 'b', E_ylim)
axes[0].set_title("E Values for Each Seed Run")

plot_with_stats(axes[1], seeds, runtime_values, "Runtime (s)", 'g', runtime_ylim)
axes[1].set_title("Runtime for Each Seed Run")

plot_with_stats(axes[2], seeds, speed_values, "Speed", 'r', speed_ylim)
axes[2].set_title("Speed for Each Seed Run")
axes[2].set_xlabel("Seed")

axes[2].set_xticks(seeds[::5])

plt.tight_layout()

plt.savefig('data/graphs.png')

plt.show()
