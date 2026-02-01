from __future__ import annotations

import math


def sat_vapor_pressure(ta_c: float) -> float:
    return 610.94 * math.exp((17.625 * ta_c) / (ta_c + 243.04))


def dew_point_from_rh(ta_c: float, rh_pct: float) -> float:
    if not math.isfinite(rh_pct) or rh_pct <= 0:
        return float("nan")
    gamma = math.log(rh_pct / 100) + (17.625 * ta_c) / (243.04 + ta_c)
    return (243.04 * gamma) / (17.625 - gamma)


def humidity_ratio(p_a: float, p_atm: float = 101325.0) -> float:
    return 0.62198 * p_a / (p_atm - p_a)
