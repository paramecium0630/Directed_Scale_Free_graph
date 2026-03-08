import numpy as np
import matplotlib.pyplot as plt

datain = np.loadtxt("output/indegree_distribution.txt", dtype=float)
dataout = np.loadtxt("output/outdegree_distribution.txt", dtype=float)

kin = datain[:, 0]
pk_in = datain[:, 1]
kout = dataout[:, 0]
pk_out = dataout[:, 1]

min_kin = min(kin)
min_kout = min(kout)
max_kin = max(kin)
max_kout = max(kout)

# 移除 k=0（因為 log(0) 不存在）
mask = (pk_in > 0) & (kin > min_kin*10) & (kin < max_kin*0.4)
kin = kin[mask]
pk_in = pk_in[mask]
mask = (pk_out > 0) & (kout > min_kout*15) & (kout < max_kout*0.8)
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
plt.figure(figsize=(8, 6))
plt.scatter(kin, pk_in, label="In-Degree Distribution", alpha=0.8)
plt.scatter(kout, pk_out, label="Out-Degree Distribution", alpha=0.8)
plt.tick_params(axis='both', which='major', labelsize=14)  # 主刻度字體
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Degree (k)", fontsize=20)
plt.ylabel("P(k)", fontsize=20)
plt.title("In-Degree and Out-Degree Distribution")

fit_pk_in = 10 ** (coeffsin[1] + coeffsin[0] * log_kin)
fit_pk_out = 10 ** (coeffsout[1] + coeffsout[0] * log_kout)
plt.plot(kin, fit_pk_in, color='tab:blue', linestyle='--', label=r"$\gamma_{{in}}$ = {:.2f}".format(-slopein))
plt.plot(kout, fit_pk_out, color='tab:orange', linestyle='--', label=r"$\gamma_{{out}}$ = {:.2f}".format(-slopeout))

plt.legend()
plt.tight_layout()
plt.savefig("figure/degree_distribution.png", dpi=300)
plt.show()
