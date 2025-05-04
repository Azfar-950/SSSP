import matplotlib.pyplot as plt

threads = [1, 2, 4, 8]
times = [10.5, 6.2, 3.8, 2.9]  # Replace with actual measured times

plt.figure(figsize=(8, 6))
plt.plot(threads, times, 'b-o', label='Execution Time')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (seconds)')
plt.title('Scalability of SSSP Update Algorithm')
plt.grid(True)
plt.legend()
plt.savefig('scalability.png')
plt.close()