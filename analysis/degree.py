import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/degree.txt", dtype=int)

N = len(data)
indeg = data[:, 0]
outdeg = data[:, 1]

distribution_in = {}
distribution_out = {}   

# Calculate degree distributions
for k in indeg:
    distribution_in[k] = distribution_in.get(k, 0) + 1
for k in outdeg:
    distribution_out[k] = distribution_out.get(k, 0) + 1

# fitting
xin = np.array(sorted(distribution_in.keys())) # 取出 distribution_in 的 key，並排序後轉成 numpy array
yin = np.array([distribution_in[k] for k in xin]) # 取出 distribution_in 裡面 key = k 的值，並轉成 numpy array
xout = np.array(sorted(distribution_out.keys())) # 取出 distribution_out 的 key，並排序後轉成 numpy array
yout = np.array([distribution_out[k] for k in xout]) # 取出 distribution_out 裡面 key = k 的值，並轉成 numpy array

# 移除 k=0（因為 log(0) 不存在）
mask = (xin > 40) & (xin < 200)
xin = xin[mask]
yin = yin[mask]
mask = (xout > 40) & (xout < 200)
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
# 對 out-degree 分布進行線性回歸
coeffsout = np.polyfit(log_xout, log_yout, 1)
slopeout = coeffsout[0]
print(f"In-degree distribution slope: {-slopein:.4f}")
print(f"Out-degree distribution slope: {-slopeout:.4f}")

# plot  
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.scatter(xin, yin, label="In-Degree Data")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("In-Degree (k)")
plt.ylabel("Count")
plt.title("In-Degree Distribution")
fit_yin = 10**(coeffsin[1] + coeffsin[0] * log_xin)
plt.plot(xin, fit_yin, color='red', label=f"Fit: slope={slopein:.2f}")
plt.legend()
plt.subplot(1, 2, 2)
plt.scatter(xout, yout, label="Out-Degree Data")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Out-Degree (k)")
plt.ylabel("Count")
plt.title("Out-Degree Distribution")
fit_yout = 10**(coeffsout[1] + coeffsout[0] * log_xout)
plt.plot(xout, fit_yout, color='red', label=f"Fit: slope={slopeout:.2f}")
plt.legend()
plt.tight_layout()
plt.show()