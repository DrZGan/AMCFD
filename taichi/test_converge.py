"""
Unit test for converge.enhance_converge_speed.

Constructs a small uniform case where the aggregated line system is
predictable, runs enhance_converge_speed, and verifies enthalpy updates
against a reference TDMA computed in Python. Uses Taichi kernels to
apply the expected correction efficiently across j,k loops.
"""

import taichi as ti
import numpy as np

from data_structures import State, DiscretCoeffs, GridParams
from converge import enhance_converge_speed


def _build_uniform_coeffs(ni, nj, nk, a=2.0, b=0.5, c=0.25, s=1.0):
    coeffs = DiscretCoeffs(ni, nj, nk)
    # Create numpy arrays
    ap = np.full((ni, nj, nk), a, dtype=np.float64)
    ae = np.full((ni, nj, nk), b, dtype=np.float64)
    aw = np.full((ni, nj, nk), c, dtype=np.float64)
    an = np.zeros((ni, nj, nk), dtype=np.float64)
    as_ = np.zeros((ni, nj, nk), dtype=np.float64)
    at = np.zeros((ni, nj, nk), dtype=np.float64)
    ab = np.zeros((ni, nj, nk), dtype=np.float64)
    su = np.full((ni, nj, nk), s, dtype=np.float64)
    # Push to fields
    coeffs.ap.from_numpy(ap)
    coeffs.ae.from_numpy(ae)
    coeffs.aw.from_numpy(aw)
    coeffs.an.from_numpy(an)
    coeffs.as_.from_numpy(as_)
    coeffs.at.from_numpy(at)
    coeffs.ab.from_numpy(ab)
    coeffs.su.from_numpy(su)
    return coeffs


def _tdma_reference(ni, nj, nk, a, b, c, s):
    """Compute reference delh[i] for uniform coefficients and zero enthalpy.
    Aggregations over j,k simply multiply by interior count.
    """
    # Interior bounds
    ist, ien = 1, ni - 2
    jst, jen = 1, nj - 2
    kst, ken = 1, nk - 2
    M = (jen - jst + 1) * (ken - kst + 1)  # interior j,k count

    bl = np.zeros(ni, dtype=np.float64)
    blp = np.zeros(ni, dtype=np.float64)
    blm = np.zeros(ni, dtype=np.float64)
    blc = np.zeros(ni, dtype=np.float64)

    # Uniform over interior i
    for i in range(ist, ien + 1):
        bl[i] = M * a
        blp[i] = M * b
        blm[i] = M * c
        blc[i] = M * s  # all h terms are zero initially

    # TDMA forward/backward
    eps = 1e-20
    pib = np.zeros(ni, dtype=np.float64)
    qib = np.zeros(ni, dtype=np.float64)
    delh = np.zeros(ni, dtype=np.float64)

    i_seed = 2
    pib[i_seed] = blp[i_seed] / (bl[i_seed] + eps)
    qib[i_seed] = blc[i_seed] / (bl[i_seed] + eps)

    for i in range(3, ni - 1):
        denom = bl[i] - blm[i] * pib[i - 1]
        denom = denom if abs(denom) > eps else eps
        pib[i] = blp[i] / denom
        qib[i] = (blc[i] + blm[i] * qib[i - 1]) / denom

    delh[ni - 2] = qib[ni - 2]
    for i in range(ni - 3, 1, -1):
        delh[i] = pib[i] * delh[i + 1] + qib[i]

    return delh


@ti.data_oriented
class _ExpectedField:
    def __init__(self, ni, nj, nk):
        self.h = ti.field(dtype=ti.f64, shape=(ni, nj, nk))
        self.delh = ti.field(dtype=ti.f64, shape=(ni,))
        self.ni = ni
        self.nj = nj
        self.nk = nk

    @ti.kernel
    def clear(self):
        for i, j, k in self.h:
            self.h[i, j, k] = 0.0

    @ti.kernel
    def apply_delh(self):
        for i, j, k in self.h:
            if 1 <= i <= self.ni - 2 and 1 <= j <= self.nj - 2 and 1 <= k <= self.nk - 2:
                self.h[i, j, k] += self.delh[i]


def test_enhance_converge_speed_uniform():
    ti.init(arch=ti.cpu)

    # Small grid with interior
    ni, nj, nk = 6, 5, 4
    a, b, c, s = 2.0, 0.5, 0.25, 1.0

    # State and coeffs
    state = State(ni, nj, nk)
    coeffs = _build_uniform_coeffs(ni, nj, nk, a=a, b=b, c=c, s=s)

    # Zero enthalpy initially
    state.enthalpy.from_numpy(np.zeros((ni, nj, nk), dtype=np.float64))

    # Run convergence acceleration
    enhance_converge_speed(state, coeffs)

    # Reference delh
    delh = _tdma_reference(ni, nj, nk, a, b, c, s)

    # Build expected enthalpy using Taichi kernel for speed
    exp = _ExpectedField(ni, nj, nk)
    exp.clear()
    exp.delh.from_numpy(delh)
    exp.apply_delh()

    h_expected = exp.h.to_numpy()
    h_result = state.enthalpy.to_numpy()

    # Compare
    assert np.allclose(h_result, h_expected, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    # Run as a simple script
    test_enhance_converge_speed_uniform()
    print("test_enhance_converge_speed_uniform: PASSED")
