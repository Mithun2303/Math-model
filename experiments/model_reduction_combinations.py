"""Phase 4, Task 3 (part 2): which parameter COMBINATIONS are recovered (reads the replicate fits of model_reduction.py; no refit).

Run after experiments/model_reduction.py. The low-burden linear block M = [[r-q-beta, a],[q, -(a+delta_D)]] has eigenvalues
lambda_+ (dominant growth) and lambda_-; m = r - q - beta is its (1,1) entry. SYNTHETIC data.
"""
from _p4common import *  # noqa: F401,F403


def lam(r, q, a, dd, b):
    M = np.array([[r - q - b, a], [q, -(a + dd)]])
    ev = np.sort(np.linalg.eigvals(M).real)[::-1]
    return ev


def main():
    d = pd.read_csv(TAB / "p4_model_reduction_replicates.csv", comment="#")
    d = d[d.model.str.startswith("M9")]
    rows = []
    for (gt, obs), g in d.groupby(["truth", "observed"]):
        th = truth(gt)
        lt = lam(th["r"], th["q"], th["a"], th["delta_D"], th["beta"])
        m_true = th["r"] - th["q"] - th["beta"]
        errs = []
        for _, r in g.iterrows():
            est = {k: th[k] * np.exp(r[f"err_{k}"]) for k in ["r", "q", "a", "delta_D", "beta"]}
            l = lam(est["r"], est["q"], est["a"], est["delta_D"], est["beta"])
            m = est["r"] - est["q"] - est["beta"]
            gg = est["r"] - est["beta"]
            errs.append({"lambda_plus": np.log(l[0] / lt[0]) if lt[0] > 0 and l[0] > 0 else np.nan,
                         "lambda_minus": np.log(abs(l[1] / lt[1])), "m=r-q-beta": np.log(m / m_true) if (m > 0 and m_true > 0) else np.nan,
                         "g=r-beta": np.log(gg / (th["r"] - th["beta"])) if gg > 0 else np.nan, "N0": r["err_N0"], "r": r["err_r"], "beta": r["err_beta"],
                         "q": r["err_q"], "a": r["err_a"], "delta_D": r["err_delta_D"]})
        e = pd.DataFrame(errs)
        row = {"truth": gt, "observed": obs, "n_rep": len(g), "lambda_plus_true": float(lt[0]), "lambda_minus_true": float(lt[1])}
        for c in e.columns:
            row[f"median_abs_ln_err_{c}"] = float(np.nanmedian(np.abs(e[c])))
        rows.append(row)
    out = pd.DataFrame(rows)
    save_table(out, "p4_model_reduction_combinations")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
