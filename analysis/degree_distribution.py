import numpy as np
import matplotlib.pyplot as plt


def weighted_linear_fit(k, pk):
    """Fit log10 P(k) = a*log10(k) + b."""
    log_k = np.log10(k)
    log_pk = np.log10(pk)
    a, b = np.polyfit(log_k, log_pk, 1)
    return a, b


def mle_exponent_from_pk(k, pk, kmin):
    """Continuous-tail MLE on binned degree PDF."""
    mask = (k >= kmin) & (pk > 0)
    if np.count_nonzero(mask) < 5:
        return np.nan, 0.0, 0
    kt = k[mask]
    wt = pk[mask]
    s = np.sum(wt * np.log(kt / (kmin - 0.5)))
    n = np.sum(wt)
    if s <= 0:
        return np.nan, n, np.count_nonzero(mask)
    tau = 1.0 + n / s
    return tau, n, np.count_nonzero(mask)


def ks_for_kmin(k, pk, kmin, tau):
    """KS distance between empirical tail CDF and fitted power-law tail CDF."""
    mask = (k >= kmin) & (pk > 0)
    kt = k[mask]
    wt = pk[mask]
    wt = wt / np.sum(wt)
    emp_cdf = np.cumsum(wt)
    model_pdf = kt ** (-tau)
    model_pdf = model_pdf / np.sum(model_pdf)
    model_cdf = np.cumsum(model_pdf)
    return np.max(np.abs(emp_cdf - model_cdf))


def auto_mle_with_kmin_scan(k, pk, kmin_start=5):
    """Scan kmin and return (best_tau, best_kmin, best_ks, tail_mass, n_bins)."""
    kmax = int(np.max(k))
    best = None
    for kmin in range(kmin_start, kmax + 1):
        tau, tail_mass, n_bins = mle_exponent_from_pk(k, pk, kmin)
        if np.isnan(tau):
            continue
        if n_bins < 20 or tail_mass < 0.01:
            continue
        ks = ks_for_kmin(k, pk, kmin, tau)
        if best is None or ks < best[2]:
            best = (tau, kmin, ks, tail_mass, n_bins)
    return best


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
mask_in_fit = (pk_in > 0) & (kin > 80) & (kin < 500)
mask_out_fit = (pk_out > 0) & (kout > 80) & (kout < 400)
kin_fit = kin[mask_in_fit]
pk_in_fit = pk_in[mask_in_fit]
kout_fit = kout[mask_out_fit]
pk_out_fit = pk_out[mask_out_fit]

# 線性擬合（舊方法）
slopein, interceptin = weighted_linear_fit(kin_fit, pk_in_fit)
slopeout, interceptout = weighted_linear_fit(kout_fit, pk_out_fit)
print(f"[log-log linear fit] In-degree exponent: {-slopein:.4f} (fit range 80<k<500)")
print(f"[log-log linear fit] Out-degree exponent: {-slopeout:.4f} (fit range 80<k<400)")

# 自動 kmin 掃描 + MLE（較穩健）
in_best = auto_mle_with_kmin_scan(kin, pk_in, kmin_start=5)
out_best = auto_mle_with_kmin_scan(kout, pk_out, kmin_start=5)
if in_best is not None:
    tau, kmin, ks, tail_mass, n_bins = in_best
    print(f"[auto MLE] In-degree exponent: {tau:.4f}, kmin={kmin}, KS={ks:.4f}, tail_mass={tail_mass:.4f}, bins={n_bins}")
else:
    print("[auto MLE] In-degree exponent: unavailable (insufficient tail data)")
if out_best is not None:
    tau, kmin, ks, tail_mass, n_bins = out_best
    print(f"[auto MLE] Out-degree exponent: {tau:.4f}, kmin={kmin}, KS={ks:.4f}, tail_mass={tail_mass:.4f}, bins={n_bins}")
else:
    print("[auto MLE] Out-degree exponent: unavailable (insufficient tail data)")

# plot
plt.figure(figsize=(8, 6))
plt.scatter(kin_fit, pk_in_fit, label="In-Degree Distribution", alpha=0.8)
plt.scatter(kout_fit, pk_out_fit, label="Out-Degree Distribution", alpha=0.8)
plt.tick_params(axis='both', which='major', labelsize=14)  # 主刻度字體
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Degree (k)", fontsize=20)
plt.ylabel("P(k)", fontsize=20)
plt.title("In-Degree and Out-Degree Distribution")

fit_pk_in = 10 ** (interceptin + slopein * np.log10(kin_fit))
fit_pk_out = 10 ** (interceptout + slopeout * np.log10(kout_fit))
plt.plot(kin_fit, fit_pk_in, color='tab:blue', linestyle='--', label=r"$\gamma_{{in,lin}}$ = {:.2f}".format(-slopein))
plt.plot(kout_fit, fit_pk_out, color='tab:orange', linestyle='--', label=r"$\gamma_{{out,lin}}$ = {:.2f}".format(-slopeout))

plt.legend()
plt.tight_layout()
plt.savefig("figure/degree_distribution.png", dpi=300)
# if "agg" not in plt.get_backend().lower():
plt.show()
