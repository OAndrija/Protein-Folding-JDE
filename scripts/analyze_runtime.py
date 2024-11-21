import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file containing the runtime results
data = pd.read_csv('data/timing_results.csv')

# Extract the threads and runtime columns
threads = data['threads']
runtime = data['time']

# Plot runtime vs. number of threads
plt.plot(threads, runtime, marker='o', linestyle='-', color='b', label='Runtime')

# Add labels and title
plt.xlabel('Number of Threads')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime vs. Number of Threads')

# Show grid for better readability
plt.grid(True)

# Save the plot as a PNG image
plt.savefig('data/runtime_analysis.png')

# Display the plot
plt.legend()
plt.show()
