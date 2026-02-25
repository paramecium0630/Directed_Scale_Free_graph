import numpy as np
import matplotlib.pyplot as plt

datain = np.loadtxt("output/indegree_distribution.txt", dtype=float)
dataout = np.loadtxt("output/outdegree_distribution.txt", dtype=float)

kin = datain[:, 0]
pk_in = datain[:, 1]
kout = dataout[:, 0]
pk_out = dataout[:, 1]

# 移除 k=0（因為 log(0) 不存在）
mask = (pk_in > 0) & (kin > 10) & (kin < 150)
kin = kin[mask]
pk_in = pk_in[mask]
mask = (pk_out > 0) & (kout > 10) & (kout < 150)
kout = kout[mask]
pk_out = pk_out[mask]

# 轉成 log
log_kin = np.log10(kin)
log_pk_in = np.log10(pk_in)
log_kout = np.log10(kout)
log_pk_out = np.log10(pk_out)

# 對 in-degree 分布進行線性回歸
coeffsin = np.polyfit(log_kin, log_pk_in, 1)
slopein = coeffsin[0]
# 對 out-degree 分布進行線性回歸
coeffsout = np.polyfit(log_kout, log_pk_out, 1)
slopeout = coeffsout[0]
print(f"In-degree distribution slope: {-slopein:.4f}")
print(f"Out-degree distribution slope: {-slopeout:.4f}")

# plot
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.scatter(kin, pk_in, label="In-Degree Distribution")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("In-Degree (k)")
plt.ylabel("P(k)")
plt.title("In-Degree Distribution")
fit_pk_in = 10**(coeffsin[1] + coeffsin[0] * log_kin)
plt.plot(kin, fit_pk_in, color='red')
plt.legend()
plt.subplot(1, 2, 2)            
plt.scatter(kout, pk_out, label="Out-Degree Distribution")
plt.xscale('log')
plt.yscale('log')   
plt.xlabel("Out-Degree (k)")
plt.ylabel("P(k)")
plt.title("Out-Degree Distribution")
fit_pk_out = 10**(coeffsout[1] + coeffsout[0] * log_kout)
plt.plot(kout, fit_pk_out, color='red')
plt.legend()
plt.tight_layout()
plt.show()