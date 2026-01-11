from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

# Proven optimal LABS energies H*(N) = min_S sum_{g=1..N-1} C_g(S)^2
# Source: OEIS A102780 / b102780.txt (proven for N <= 66).
# Note: N starts at 1, and H*(1)=0.
_OEIS_A102780_HSTAR = {
    1: 0,  2: 1,  3: 1,  4: 2,  5: 2,  6: 7,  7: 3,  8: 8,  9: 12, 10: 13,
    11: 5, 12: 10, 13: 6, 14: 19, 15: 15, 16: 24, 17: 32, 18: 25, 19: 29, 20: 26,
    21: 26, 22: 39, 23: 47, 24: 36, 25: 36, 26: 45, 27: 37, 28: 50, 29: 62, 30: 59,
    31: 67, 32: 64, 33: 64, 34: 65, 35: 73, 36: 82, 37: 86, 38: 87, 39: 99, 40: 108,
    41: 108, 42: 101, 43: 109, 44: 122, 45: 118, 46: 131, 47: 135, 48: 140, 49: 136, 50: 153,
    51: 153, 52: 166, 53: 170, 54: 175, 55: 171, 56: 192, 57: 188, 58: 197, 59: 205, 60: 218,
    61: 226, 62: 235, 63: 207, 64: 208, 65: 240, 66: 257,
}

def labs_energy(sequence_pm1: List[int]) -> int:
    """
    LABS / Bernasconi open-boundary energy:
        H(S) = sum_{g=1..N-1} C_g^2
        C_g = sum_{i=0..N-g-1} s_i s_{i+g},  s_i in {+1,-1}
    """
    s = sequence_pm1
    n = len(s)
    e = 0
    for g in range(1, n):
        cg = 0
        for i in range(n - g):
            cg += s[i] * s[i + g]
        e += cg * cg
    return e

def merit_factor(n: int, energy: int) -> float:
    """
    Merit factor F = n^2 / (2 * H). (Undefined/infinite for H=0.)
    """
    if energy == 0:
        return float("inf")
    return (n * n) / (2.0 * energy)

@dataclass
class LabsOptimumResult:
    n: int
    energy: int
    proven: bool                    # True if from the proven table or solver proved optimality
    status: str                     # "KNOWN_PROVEN", "OPTIMAL", "FEASIBLE", "UNKNOWN", ...
    sequence_pm1: Optional[List[int]] = None  # Returned only if solved via CP-SAT

def labs_true_optimal_energy(
    n: int,
    *,
    method: str = "auto",           # "auto" | "lookup" | "solve"
    time_limit_s: Optional[float] = 60.0,
    num_workers: int = 8,
    symmetry_breaking: bool = True,
) -> LabsOptimumResult:
    """
    Returns the proven optimal LABS energy if known (N<=66).
    Otherwise attempts an exact solve with OR-Tools CP-SAT.

    If CP-SAT returns status OPTIMAL, the result is "true optimal".
    If FEASIBLE/UNKNOWN, it is NOT guaranteed optimal.
    """
    if n <= 0:
        raise ValueError("N must be a positive integer.")
    if method not in {"auto", "lookup", "solve"}:
        raise ValueError("method must be one of: 'auto', 'lookup', 'solve'.")

    # 1) Proven lookup (fast, exact)
    if n in _OEIS_A102780_HSTAR:
        return LabsOptimumResult(
            n=n,
            energy=_OEIS_A102780_HSTAR[n],
            proven=True,
            status="KNOWN_PROVEN",
            sequence_pm1=None,
        )

    if method == "lookup":
        raise ValueError(
            f"No proven optimal value is available in the built-in table for N={n} "
            f"(table covers N<=66). Use method='solve' to attempt an exact solve."
        )

    # 2) Attempt exact solve with CP-SAT
    try:
        from ortools.sat.python import cp_model
    except Exception as exc:
        raise RuntimeError(
            "OR-Tools is required for method='solve'. Install with:\n"
            "  pip install ortools\n"
        ) from exc

    model = cp_model.CpModel()

    # Bits x_i in {0,1} representing spins s_i in {+1,-1} via s_i = 1 - 2*x_i
    x = [model.NewBoolVar(f"x_{i}") for i in range(n)]

    # Optional symmetry breaking: fix first spin to +1 (x_0 = 0) to remove global complement symmetry
    if symmetry_breaking:
        model.Add(x[0] == 0)

    # Build C_g and C_g^2 for each g
    c_vars = []
    e_vars = []
    for g in range(1, n):
        m = n - g  # number of terms in C_g
        # eq_i = 1 iff x_i == x_{i+g}
        eq = [model.NewBoolVar(f"eq_{i}_{g}") for i in range(m)]
        for i in range(m):
            model.Add(x[i] == x[i + g]).OnlyEnforceIf(eq[i])
            model.Add(x[i] != x[i + g]).OnlyEnforceIf(eq[i].Not())

        # C_g = sum_i ( +1 if equal else -1 ) = sum_i (2*eq_i - 1)
        c_min, c_max = -m, m
        Cg = model.NewIntVar(c_min, c_max, f"C_{g}")
        model.Add(Cg == sum((2 * eq_i - 1) for eq_i in eq))
        c_vars.append(Cg)

        # E_g = C_g^2
        Eg = model.NewIntVar(0, m * m, f"E_{g}")
        model.AddMultiplicationEquality(Eg, [Cg, Cg])
        e_vars.append(Eg)

    H = model.NewIntVar(0, sum((n - g) * (n - g) for g in range(1, n)), "H")
    model.Add(H == sum(e_vars))
    model.Minimize(H)

    solver = cp_model.CpSolver()
    if time_limit_s is not None:
        solver.parameters.max_time_in_seconds = float(time_limit_s)
    solver.parameters.num_search_workers = int(num_workers)

    status = solver.Solve(model)
    status_name = solver.StatusName(status)

    if status_name not in {"OPTIMAL", "FEASIBLE"}:
        return LabsOptimumResult(n=n, energy=-1, proven=False, status=status_name, sequence_pm1=None)

    energy = int(solver.Value(H))
    seq_pm1 = [1 - 2 * int(solver.Value(xi)) for xi in x]  # map bits -> +/-1

    # If OPTIMAL, this is a proven optimum for that N (within the model definition)
    proven = (status_name == "OPTIMAL")
    return LabsOptimumResult(n=n, energy=energy, proven=proven, status=status_name, sequence_pm1=seq_pm1)


if __name__ == "__main__":
    # Examples
    for N in [30, 66]:
        r = labs_true_optimal_energy(N, method="auto")
        print(f"N={N}: H*={r.energy}, proven={r.proven}, status={r.status}, F*={merit_factor(N, r.energy)}")

    # For N>66: attempt an exact solve (may be slow!)
    r = labs_true_optimal_energy(70, method="solve", time_limit_s=30.0, num_workers=8)
    print(f"N=70: H={r.energy}, proven={r.proven}, status={r.status}")
    if r.sequence_pm1:
        print("Sequence (±1):", r.sequence_pm1)
        print("Check energy:", labs_energy(r.sequence_pm1))

#dummy
