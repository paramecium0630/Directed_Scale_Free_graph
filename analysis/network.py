import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

data_adj = np.loadtxt("output/directed_scale_free_graph.txt", dtype=int)

N = int(np.max(data_adj))

# Reshape Aij into an adjacency matrix
Aij = np.zeros((N, N))
for row in data_adj:
    i, j = row[0]-1, row[1]-1  # Convert to 0-based indexing
    Aij[i, j] = 1

Aij = Aij.T # Transpose to match the convention A(i,j) = 1 if j->i

# Create a directed graph from the adjacency matrix
G = nx.from_numpy_array(Aij, create_using=nx.DiGraph)

# Calculate in-degrees and out-degrees
in_degrees = dict(G.in_degree())
out_degrees = dict(G.out_degree())

distribution_in = {}
distribution_out = {}

# Calculate degree distributions
for k in in_degrees.values():
    # 從字典 distribution_in 裡面「取出 key = k 的值」
    distribution_in[k] = distribution_in.get(k, 0) + 1 
for k in out_degrees.values():
    distribution_out[k] = distribution_out.get(k, 0) + 1

# fitting

xin = np.array(sorted(distribution_in.keys())) # 取出 distribution_in 的 key，並排序後轉成 numpy array
yin = np.array([distribution_in[k] for k in xin]) # 取出 distribution_in 裡面 key = k 的值，並轉成 numpy array
xout = np.array(sorted(distribution_out.keys())) # 取出 distribution_out 的 key，並排序後轉成 numpy array
yout = np.array([distribution_out[k] for k in xout]) # 取出 distribution_out 裡面 key = k 的值，並轉成 numpy array

# 移除 k=0（因為 log(0) 不存在）

xmin = 90; xmax = 300

mask = (xin > xmin) & (xin < xmax)
xin = xin[mask]
yin = yin[mask]
mask = (xout > xmin) & (xout < xmax)
xout = xout[mask]
yout = yout[mask]

# 轉成 log
log_xin = np.log10(xin)
log_yin = np.log10(yin)
log_xout = np.log10(xout)
log_yout = np.log10(yout)

# 對 in-degree 分布進行線性回歸
coeffsin = np.polyfit(log_xin, log_yin, 1)
slopein = coeffsin[0]
interceptin = coeffsin[1]
gammain = -slopein
print("Estimated gamma =", gammain)

# 對 out-degree 分布進行線性回歸
coeffsout = np.polyfit(log_xout, log_yout, 1)
slopeout = coeffsout[0]
interceptout = coeffsout[1]
gammaout = -slopeout
print("Estimated gamma =", gammaout)

# Plot in-degree distribution
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.scatter(xin, yin, label="In-Degree Data")
plt.xscale('log')
plt.yscale('log')

# 畫回 power-law 線
fit_yin = 10**(interceptin + slopein * log_xin)
plt.plot(xin, fit_yin, color='red', label=f"Fit γ={slopein:.2f}")

plt.xlabel('In-Degree (k_in)')
plt.ylabel('Number of Nodes')
plt.title('In-Degree Distribution')

# Plot out-degree distribution
plt.subplot(1, 2, 2)
plt.scatter(xout, yout, label="Out-Degree Data", color='orange')
plt.xscale('log')
plt.yscale('log')

# 畫回 power-law 線
fit_yout = 10**(interceptout + slopeout * log_xout)
plt.plot(xout, fit_yout, color='red', label=f"Fit γ={slopeout:.2f}")

plt.xlabel('Out-Degree (k_out)')
plt.ylabel('Number of Nodes')
plt.title('Out-Degree Distribution')

plt.savefig("figure/degree_distributions.png")
plt.show()