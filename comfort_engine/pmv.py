from __future__ import annotations

import math
from typing import Dict, List

from .constants import CLO_TO_ICL, EPSILON, F_EFF, MET_TO_WM2, SIGMA
from .psychrometrics import sat_vapor_pressure
from .types import ComfortState, ModelConfig, Occupant, PMVResult, TraceStep


def _is_finite(value: float) -> bool:
    return math.isfinite(value)


def compute_pmv(
    state: ComfortState,
    occupant: Occupant,
    config: ModelConfig | None = None,
) -> PMVResult:
    cfg = config or ModelConfig()
    trace: List[TraceStep] = []
    warnings: List[str] = []

    def add_trace(
        step_id: str,
        equation: str,
        inputs: Dict[str, float] | None = None,
        outputs: Dict[str, float] | None = None,
        units: Dict[str, str] | None = None,
        condition: str | None = None,
    ) -> None:
        trace.append(
            TraceStep(
                id=step_id,
                equation=equation,
                inputs=inputs or {},
                outputs=outputs or {},
                units=units or {},
                condition=condition,
            )
        )

    Ta = state.Ta
    Tr = state.Tr
    RH = state.RH
    v = state.v if _is_finite(state.v) else 0.0
    patm = state.patm

    m_wm2 = occupant.m_wm2
    clo = occupant.clo
    wme = occupant.wme

    if not _is_finite(m_wm2):
        raise ValueError("Occupant.m_wm2 must be finite.")
    if not _is_finite(clo):
        raise ValueError("Occupant.clo must be finite.")

    M = m_wm2
    H = M - wme
    met_equiv = M / MET_TO_WM2 if MET_TO_WM2 else float("nan")
    add_trace(
        "met_input",
        "M provided in W/m2; met_equiv = M / 58.15",
        inputs={"M": M},
        outputs={"H": H, "met_equiv": met_equiv},
        units={"M": "W/m2", "H": "W/m2", "met_equiv": "met"},
    )

    Icl = clo * CLO_TO_ICL
    f_cl = 1 + 1.29 * Icl if Icl <= 0.078 else 1.05 + 0.645 * Icl
    add_trace(
        "clothing",
        "Icl = clo * 0.155; f_cl = 1 + 1.29*Icl (<=0.078) else 1.05 + 0.645*Icl",
        inputs={"clo": clo},
        outputs={"Icl": Icl, "f_cl": f_cl},
        units={"Icl": "m2K/W"},
        condition="Icl <= 0.078" if Icl <= 0.078 else "Icl > 0.078",
    )

    p_ws = sat_vapor_pressure(Ta)
    p_a = (RH / 100.0) * p_ws
    add_trace(
        "humidity",
        "p_ws = 610.94 * exp(17.625*Ta/(Ta+243.04)); p_a = RH/100 * p_ws",
        inputs={"Ta": Ta, "RH": RH},
        outputs={"p_ws": p_ws, "p_a": p_a},
        units={"p_ws": "Pa", "p_a": "Pa"},
    )

    h_forced = 12.1 * math.sqrt(max(v, 0.0))
    taa = Ta + 273.0
    tra = Tr + 273.0
    tcla = taa + (35.5 - Ta) / (3.5 * Icl + 0.1)

    p1 = Icl * f_cl
    p2 = p1 * 3.96
    p3 = p1 * 100.0
    p4 = p1 * taa
    p5 = 308.7 - 0.028 * H + p2 * (tra / 100.0) ** 4

    add_trace(
        "tcl_seed",
        "tcla = Ta+273 + (35.5-Ta)/(3.5*Icl+0.1)",
        inputs={"Ta": Ta, "Icl": Icl},
        outputs={"tcla": tcla},
        units={"tcla": "K"},
    )

    xn = tcla / 100.0
    xf = xn
    h_c = h_forced
    iterations = 0
    residual = float("inf")

    if cfg.convection not in {"auto", "forced", "natural"}:
        raise ValueError("ModelConfig.convection must be auto, forced, or natural.")

    for i in range(cfg.max_iter):
        xf = (xf + xn) / 2.0
        h_natural_iter = 2.38 * math.pow(abs(100.0 * xf - taa), 0.25)
        if cfg.convection == "auto":
            h_c = max(h_forced, h_natural_iter)
        elif cfg.convection == "forced":
            h_c = h_forced
        else:
            h_c = h_natural_iter
        xn = (p5 + p4 * h_c - p2 * (xf ** 4)) / (100.0 + p3 * h_c)
        residual = abs(xn - xf)
        iterations = i + 1
        if residual <= cfg.tolerance:
            break

    Tcl = 100.0 * xn - 273.0
    h_natural = 2.38 * math.pow(abs(Tcl - Ta), 0.25)
    if cfg.convection == "auto":
        h_c = max(h_forced, h_natural)
    elif cfg.convection == "forced":
        h_c = h_forced
    else:
        h_c = h_natural

    add_trace(
        "tcl_solver",
        "Iterate Tcl until |xn-xf| <= tolerance",
        inputs={"max_iter": cfg.max_iter, "tolerance": cfg.tolerance},
        outputs={
            "Tcl": Tcl,
            "h_forced": h_forced,
            "h_natural": h_natural,
            "h_c": h_c,
            "iterations": float(iterations),
            "residual": residual,
        },
        units={"Tcl": "C", "h_c": "W/m2K"},
    )

    Tcl_K = Tcl + 273.15
    Tr_K = Tr + 273.15
    q_rad = 3.96e-8 * f_cl * (Tcl_K ** 4 - Tr_K ** 4)
    q_conv = f_cl * h_c * (Tcl - Ta)
    add_trace(
        "dry_heat",
        "Q_rad = 3.96e-8*f_cl*(Tcl_K^4-Tr_K^4); Q_conv = f_cl*h_c*(Tcl-Ta)",
        inputs={"Tcl": Tcl, "Tr": Tr, "Ta": Ta, "f_cl": f_cl, "h_c": h_c},
        outputs={"Q_rad": q_rad, "Q_conv": q_conv},
        units={"Q_rad": "W/m2", "Q_conv": "W/m2"},
    )

    e_diff = 3.05e-3 * (5733.0 - 6.99 * H - p_a)
    e_sw = max(0.0, 0.42 * (H - MET_TO_WM2))
    e_res = 1.7e-5 * M * (5867.0 - p_a)
    c_res = 0.0014 * M * (34.0 - Ta)

    e_sk = e_diff + e_sw
    q_res = e_res + c_res
    L = q_conv + q_rad + e_sk + q_res

    add_trace(
        "latent_heat",
        "E_diff + E_sw, E_res + C_res",
        inputs={"H": H, "M": M, "p_a": p_a, "Ta": Ta},
        outputs={
            "E_diff": e_diff,
            "E_sw": e_sw,
            "E_res": e_res,
            "C_res": c_res,
            "E_sk": e_sk,
            "Q_res": q_res,
            "L": L,
        },
        units={"E_sk": "W/m2", "Q_res": "W/m2", "L": "W/m2"},
    )

    pmv_factor = 0.303 * math.exp(-0.036 * M) + 0.028
    PMV = pmv_factor * (H - L)
    PPD = 100.0 - 95.0 * math.exp(-0.03353 * (PMV ** 4) - 0.2179 * (PMV ** 2))
    add_trace(
        "pmv",
        "PMV = (0.303*exp(-0.036*M)+0.028)*(H-L); PPD = 100 - 95*exp(...) ",
        inputs={"M": M, "H": H, "L": L},
        outputs={"PMV": PMV, "PPD": PPD, "pmv_factor": pmv_factor},
    )

    delta_tr = Tcl - Tr
    hr = float("nan")
    hr_method = "direct"
    if cfg.hr_method not in {"auto", "direct", "linearized"}:
        raise ValueError("ModelConfig.hr_method must be auto, direct, or linearized.")

    if cfg.hr_method == "direct":
        if abs(delta_tr) <= 0.05:
            warnings.append("Small Tcl-Tr; direct hr may be noisy.")
        hr = q_rad / (f_cl * delta_tr) if abs(delta_tr) > 0 else float("nan")
        hr_method = "direct"
    elif cfg.hr_method == "linearized":
        t_mean = (Tcl + Tr) / 2.0 + 273.15
        hr = 4.0 * EPSILON * SIGMA * F_EFF * (t_mean ** 3)
        hr_method = "linearized"
    else:
        if _is_finite(delta_tr) and abs(delta_tr) > 0.05:
            hr = q_rad / (f_cl * delta_tr)
            hr_method = "direct"
        else:
            t_mean = (Tcl + Tr) / 2.0 + 273.15
            hr = 4.0 * EPSILON * SIGMA * F_EFF * (t_mean ** 3)
            hr_method = "linearized"

    eq_slope = -hr / h_c if _is_finite(hr) and _is_finite(h_c) else float("nan")
    dominance = "-"
    if _is_finite(hr) and _is_finite(h_c) and h_c > 0:
        ratio = hr / h_c
        if ratio > 1.2:
            dominance = "Radiation-dominant"
        elif ratio < 0.8:
            dominance = "Convection-dominant"
        else:
            dominance = "Mixed"

    add_trace(
        "radiation_linearization",
        "hr = Q_rad/(f_cl*(Tcl-Tr)) or 4*eps*sigma*F_eff*Tmean^3",
        inputs={"Tcl": Tcl, "Tr": Tr, "Q_rad": q_rad, "f_cl": f_cl},
        outputs={"hr": hr, "eq_slope": eq_slope, "dominance": dominance},
        condition=hr_method,
    )

    if v > 0.2:
        warnings.append("Elevated air speed: PMV may be unreliable; consider SET.")
    if met_equiv < 0.8 or met_equiv > 2.0:
        warnings.append("Met outside typical PMV range (0.8-2.0).")
    if clo > 2.0:
        warnings.append("High clothing insulation (>2 clo).")
    if Ta < 10 or Ta > 35 or Tr < 10 or Tr > 35:
        warnings.append("Air/radiant temperature outside typical PMV bounds.")
    if RH < 20 or RH > 80:
        warnings.append("Extreme humidity may reduce PMV reliability.")
    if residual > cfg.tolerance:
        warnings.append(f"Tcl solver not converged (residual {residual:.5f}).")

    branches = {
        "convection_regime": "forced" if h_forced >= h_natural else "natural",
        "hr_method": hr_method,
    }

    components = {
        "Q_conv": q_conv,
        "Q_rad": q_rad,
        "E_diff": e_diff,
        "E_sw": e_sw,
        "E_res": e_res,
        "C_res": c_res,
        "E_sk": e_sk,
        "Q_res": q_res,
        "L_total": L,
    }

    components_w = None
    if occupant.body_area_m2:
        components_w = {key: value * occupant.body_area_m2 for key, value in components.items()}

    solver = {
        "iterations": float(iterations),
        "residual": residual,
        "converged": residual <= cfg.tolerance,
    }

    inputs = {
        "Ta": Ta,
        "Tr": Tr,
        "RH": RH,
        "v": v,
        "patm": patm,
        "M_wm2": M,
        "met_equiv": met_equiv,
        "clo": clo,
        "wme": wme,
    }

    return PMVResult(
        PMV=PMV,
        PPD=PPD,
        Tcl=Tcl,
        h_c=h_c,
        q_conv=q_conv,
        q_rad=q_rad,
        e_sk=e_sk,
        q_res=q_res,
        L=L,
        M=M,
        H=H,
        pmv_factor=pmv_factor,
        eq_slope=eq_slope,
        dominance=dominance,
        trace=trace,
        warnings=warnings,
        branches=branches,
        solver=solver,
        components=components,
        components_w=components_w,
        inputs=inputs,
    )
