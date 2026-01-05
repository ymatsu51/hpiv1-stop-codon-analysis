# =========================================================
# Single-run forward simulation (2000 generations)
# One seed per run; repeat externally for 10 independent runs
# =========================================================

import random

STOP_CODONS = {"TAA", "TAG", "TGA"}


def count_stops(seq: str, frame: int):
    s = 0
    c = 0
    for i in range(frame, len(seq) - 2, 3):
        c += 1
        if seq[i : i + 3] in STOP_CODONS:
            s += 1
    return s, c


def mutate(seq: str, mu: float):
    bases = ["A", "C", "G", "T"]
    out = list(seq)
    for i, b in enumerate(out):
        if random.random() < mu:
            out[i] = random.choice([x for x in bases if x != b])
    return "".join(out)


def run_experiment(L=759, N=300, gens=2000, mu=2e-4, seed=1):
    L = int(L)
    N = int(N)
    gens = int(gens)
    seed = int(seed)
    mu = float(mu)

    random.seed(seed)
    bases = ["A", "C", "G", "T"]

    # random ancestor, constrained only by no STOP in main frame (frame 0)
    while True:
        anc = "".join(random.choice(bases) for _ in range(L))
        if count_stops(anc, 0)[0] == 0:
            break

    pop_sel = [anc] * N
    pop_neu = [anc] * N

    def pooled_stop(pop, frame):
        ts = 0
        tc = 0
        for x in pop:
            s_, c_ = count_stops(x, frame)
            ts += s_
            tc += c_
        return ts / tc

    for _g in range(1, gens + 1):
        # selection on main ORF only
        muts = [mutate(x, mu) for x in pop_sel]
        kept = [x for x in muts if count_stops(x, 0)[0] == 0]
        if len(kept) == 0:
            kept = pop_sel
        pop_sel = random.choices(kept, k=N)

        # neutral evolution
        muts2 = [mutate(x, mu) for x in pop_neu]
        pop_neu = random.choices(muts2, k=N)

    # pseudo-V corresponds to +1 frame
    return pooled_stop(pop_sel, 1), pooled_stop(pop_neu, 1)


def main():
    # ---- parameters (fixed) ----
    L = 759
    N = 300
    mu = 2e-4
    gens = 2000

    # ---- choose ONE seed per run ----
    seed = 1  # change this to 1, 2, 3, ..., 10 and run each simulation separately

    stop_sel, stop_neu = run_experiment(L=L, N=N, gens=gens, mu=mu, seed=seed)

    out = {
        "seed": seed,
        "stop_sel": float(stop_sel),
        "stop_neu": float(stop_neu),
        "diff": float(stop_sel - stop_neu),
        "ratio": float(stop_sel / stop_neu),
    }

    print(out)


if __name__ == "__main__":
    main()
