import os
from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import yaml
import taichi as ti

from data_structures import GridParams

ti.init(arch=ti.cpu)

def load_yaml_config(yaml_path: str) -> Dict:
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"Config not found: {yaml_path}")
    with open(yaml_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _ensure_list(value: Union[List, Tuple, int, float]) -> List[float]:
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    return [float(value)]


def _axis_params(axis_cfg: Dict) -> Tuple[List[float], List[int], List[float]]:
    zones = int(axis_cfg.get("zones", 1))
    lengths = _ensure_list(axis_cfg.get("zone_length_m", 1.0))
    counts_raw = axis_cfg.get("cv_per_zone", 10)
    exps_raw = axis_cfg.get("cv_boundary_exponent", 1)

    if len(lengths) != zones:
        if len(lengths) == 1:
            lengths = lengths * zones
        else:
            raise ValueError("zone_length_m length must match zones")

    if isinstance(counts_raw, (list, tuple)):
        counts = [int(v) for v in counts_raw]
    else:
        counts = [int(counts_raw)] * zones

    if isinstance(exps_raw, (list, tuple)):
        exps = [float(v) for v in exps_raw]
    else:
        exps = [float(exps_raw)] * zones

    if len(counts) != zones or len(exps) != zones:
        raise ValueError("cv_per_zone and cv_boundary_exponent must match zones")

    return lengths, counts, exps


def build_gridparams(cfg: Dict) -> GridParams:
    geom_cfg = cfg["geometry"]

    xzone, ncvx, powrx = _axis_params(geom_cfg["x"])  # lengths, counts, exponents
    yzone, ncvy, powry = _axis_params(geom_cfg["y"])  # lengths, counts, exponents
    zzone, ncvz, powrz = _axis_params(geom_cfg["z"])  # lengths, counts, exponents

    # X faces and centers
    ni = 2 + int(sum(ncvx))
    xu = np.zeros(ni, dtype=np.float64)
    x = np.zeros(ni, dtype=np.float64)
    ist = 2
    statloc = 0.0
    for L, N, p in zip(xzone, ncvx, powrx):
        for j in range(1, N + 1):
            r = j / float(N)
            term = (r ** p) if p >= 0.0 else 1.0 - (1.0 - r) ** (-p)
            xu[ist + j - 1] = statloc + L * term
        ist += N
        statloc += L
    nim1 = ni - 1
    for i in range(0, nim1):
        x[i] = 0.5 * (xu[i + 1] + xu[i])
    x[nim1] = xu[nim1]

    # Y faces and centers
    nj = 2 + int(sum(ncvy))
    yv = np.zeros(nj, dtype=np.float64)
    y = np.zeros(nj, dtype=np.float64)
    ist = 2
    statloc = 0.0
    for L, N, p in zip(yzone, ncvy, powry):
        for j in range(1, N + 1):
            r = j / float(N)
            term = (r ** p) if p >= 0.0 else 1.0 - (1.0 - r) ** (-p)
            yv[ist + j - 1] = statloc + L * term
        ist += N
        statloc += L
    njm1 = nj - 1
    for j in range(0, njm1):
        y[j] = 0.5 * (yv[j + 1] + yv[j])
    y[njm1] = yv[njm1]

    # Z faces and centers
    nk = 2 + int(sum(ncvz))
    zw = np.zeros(nk, dtype=np.float64)
    z = np.zeros(nk, dtype=np.float64)
    ist = 2
    statloc = 0.0
    for L, N, p in zip(zzone, ncvz, powrz):
        for j in range(1, N + 1):
            r = j / float(N)
            term = (r ** p) if p >= 0.0 else 1.0 - (1.0 - r) ** (-p)
            zw[ist + j - 1] = statloc + L * term
        ist += N
        statloc += L
    nkm1 = nk - 1
    for k in range(0, nkm1):
        z[k] = 0.5 * (zw[k + 1] + zw[k])
    z[nkm1] = zw[nkm1]

    # Derive cell widths (use face spacings; replicate boundary for last index)
    dx = np.zeros(ni, dtype=np.float64)
    dy = np.zeros(nj, dtype=np.float64)
    dz = np.zeros(nk, dtype=np.float64)
    dx[:-1] = xu[1:] - xu[:-1]
    dx[-1] = dx[-2] if ni > 1 else 0.0
    dy[:-1] = yv[1:] - yv[:-1]
    dy[-1] = dy[-2] if nj > 1 else 0.0
    dz[:-1] = zw[1:] - zw[:-1]
    dz[-1] = dz[-2] if nk > 1 else 0.0

    # Inverse distances (both west/east, south/north, bottom/top)
    dxpwinv = np.zeros(ni, dtype=np.float64)
    dxpeinv = np.zeros(ni, dtype=np.float64)
    dypsinv = np.zeros(nj, dtype=np.float64)
    dypninv = np.zeros(nj, dtype=np.float64)
    dzpbinv = np.zeros(nk, dtype=np.float64)
    dzptinv = np.zeros(nk, dtype=np.float64)

    for i in range(1, ni):
        dxpwinv[i] = 1.0 / (x[i] - x[i - 1]) if x[i] != x[i - 1] else 0.0
    for i in range(0, ni - 1):
        dxpeinv[i] = 1.0 / (x[i + 1] - x[i]) if x[i + 1] != x[i] else 0.0
    dxpeinv[-1] = dxpeinv[-2] if ni > 1 else 0.0

    for j in range(1, nj):
        dypsinv[j] = 1.0 / (y[j] - y[j - 1]) if y[j] != y[j - 1] else 0.0
    for j in range(0, nj - 1):
        dypninv[j] = 1.0 / (y[j + 1] - y[j]) if y[j + 1] != y[j] else 0.0
    dypninv[-1] = dypninv[-2] if nj > 1 else 0.0

    for k in range(1, nk):
        dzpbinv[k] = 1.0 / (z[k] - z[k - 1]) if z[k] != z[k - 1] else 0.0
    for k in range(0, nk - 1):
        dzptinv[k] = 1.0 / (z[k + 1] - z[k]) if z[k + 1] != z[k] else 0.0
    dzptinv[-1] = dzptinv[-2] if nk > 1 else 0.0

    # Interpolation fractions
    fracx = np.zeros(ni, dtype=np.float64)
    fracy = np.zeros(nj, dtype=np.float64)
    fracz = np.zeros(nk, dtype=np.float64)
    for i in range(0, ni - 1):
        denom = x[i + 1] - x[i]
        fracx[i] = (x[i + 1] - xu[i + 1]) / denom if denom != 0.0 else 0.0
    for j in range(0, nj - 1):
        denom = y[j + 1] - y[j]
        fracy[j] = (y[j + 1] - yv[j + 1]) / denom if denom != 0.0 else 0.0
    for k in range(0, nk - 1):
        denom = z[k + 1] - z[k]
        fracz[k] = (z[k + 1] - zw[k + 1]) / denom if denom != 0.0 else 0.0

    # Volumes and areas (interior nodes only, others remain zero)
    vol = np.zeros((ni, nj, nk), dtype=np.float64)
    volume_u = np.zeros((ni, nj, nk), dtype=np.float64)
    volume_v = np.zeros((ni, nj, nk), dtype=np.float64)
    volume_w = np.zeros((ni, nj, nk), dtype=np.float64)
    areaij = np.zeros((ni, nj), dtype=np.float64)
    areaik = np.zeros((ni, nk), dtype=np.float64)
    areajk = np.zeros((nj, nk), dtype=np.float64)
    areauij = np.zeros((ni, nj), dtype=np.float64)
    areauik = np.zeros((ni, nk), dtype=np.float64)
    areavij = np.zeros((ni, nj), dtype=np.float64)
    areavjk = np.zeros((nj, nk), dtype=np.float64)
    areawik = np.zeros((ni, nk), dtype=np.float64)
    areawjk = np.zeros((nj, nk), dtype=np.float64)

    nim1 = ni - 1
    njm1 = nj - 1
    nkm1 = nk - 1
    for k in range(1, nkm1):
        for j in range(1, njm1):
            for i in range(1, nim1):
                vol[i, j, k] = (xu[i + 1] - xu[i]) * (yv[j + 1] - yv[j]) * (zw[k + 1] - zw[k])
                volume_u[i, j, k] = (x[i] - x[i - 1]) * (yv[j + 1] - yv[j]) * (zw[k + 1] - zw[k])
                volume_v[i, j, k] = (xu[i + 1] - xu[i]) * (y[j] - y[j - 1]) * (zw[k + 1] - zw[k])
                volume_w[i, j, k] = (xu[i + 1] - xu[i]) * (yv[j + 1] - yv[j]) * (z[k] - z[k - 1])

    for j in range(1, njm1):
        for i in range(1, nim1):
            areaij[i, j] = (xu[i + 1] - xu[i]) * (yv[j + 1] - yv[j])
            areauij[i, j] = (x[i] - x[i - 1]) * (yv[j + 1] - yv[j])
            areavij[i, j] = (xu[i + 1] - xu[i]) * (y[j] - y[j - 1])

    for k in range(1, nkm1):
        for i in range(1, nim1):
            areaik[i, k] = (xu[i + 1] - xu[i]) * (zw[k + 1] - zw[k])
            areawik[i, k] = (xu[i + 1] - xu[i]) * (z[k] - z[k - 1])
            areauik[i, k] = (x[i] - x[i - 1]) * (zw[k + 1] - zw[k])

    for k in range(1, nkm1):
        for j in range(1, njm1):
            areajk[j, k] = (yv[j + 1] - yv[j]) * (zw[k + 1] - zw[k])
            areavjk[j, k] = (y[j] - y[j - 1]) * (zw[k + 1] - zw[k])
            areawjk[j, k] = (yv[j + 1] - yv[j]) * (z[k] - z[k - 1])

    # Create GridParams and populate (Taichi should already be initialized)
    gp = GridParams(ni=ni, nj=nj, nk=nk)

    gp.x.from_numpy(x)
    gp.y.from_numpy(y)
    gp.z.from_numpy(z)
    gp.xu.from_numpy(xu)
    gp.yv.from_numpy(yv)
    gp.zw.from_numpy(zw)

    gp.dx.from_numpy(dx)
    gp.dy.from_numpy(dy)
    gp.dz.from_numpy(dz)

    gp.dxpwinv.from_numpy(dxpwinv)
    gp.dxpeinv.from_numpy(dxpeinv)
    gp.dypsinv.from_numpy(dypsinv)
    gp.dypninv.from_numpy(dypninv)
    gp.dzpbinv.from_numpy(dzpbinv)
    gp.dzptinv.from_numpy(dzptinv)

    gp.vol.from_numpy(vol)
    gp.volume_u.from_numpy(volume_u)
    gp.volume_v.from_numpy(volume_v)
    gp.volume_w.from_numpy(volume_w)
    gp.areaij.from_numpy(areaij)
    gp.areaik.from_numpy(areaik)
    gp.areajk.from_numpy(areajk)
    gp.areauij.from_numpy(areauij)
    gp.areauik.from_numpy(areauik)
    gp.areavij.from_numpy(areavij)
    gp.areavjk.from_numpy(areavjk)
    gp.areawik.from_numpy(areawik)
    gp.areawjk.from_numpy(areawjk)
    gp.fracx.from_numpy(fracx)
    gp.fracy.from_numpy(fracy)
    gp.fracz.from_numpy(fracz)

    return gp


def build_from_yaml(yaml_path: str) -> GridParams:
    cfg = load_yaml_config(yaml_path)
    return build_gridparams(cfg)


def default_yaml_path() -> str:
    base = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(base)
    return os.path.join(root, "inputfile", "input_param.yaml")


def get_gridparams(yaml_path: Optional[str] = None) -> GridParams:
    """Convenience function to construct and return GridParams.

    If `yaml_path` is None, uses the default repo-root inputfile path.
    """
    path = yaml_path or default_yaml_path()
    return build_from_yaml(path)


__all__ = [
    "GridParams",
    "build_gridparams",
    "build_from_yaml",
    "get_gridparams",
]


if __name__ == "__main__":
    ypath = os.environ.get("AMC_YAML", default_yaml_path())
    gp = get_gridparams(ypath)
    ni, nj, nk = gp.ni, gp.nj, gp.nk
    # Print basic confirmation similar to Fortran dims
    print(f"ni={ni}, nj={nj}, nk={nk}")
    print(f"dimx={gp.x[ni-1]:.6e}, dimy={gp.y[nj-1]:.6e}, dimz={gp.z[nk-1]:.6e}")
