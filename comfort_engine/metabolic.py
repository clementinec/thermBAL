from __future__ import annotations

from typing import Dict

from .constants import MET_TO_WM2


def met_to_wm2(met: float, met_wm2: float = MET_TO_WM2) -> float:
    return met * met_wm2


def wm2_to_met(m_wm2: float, met_wm2: float = MET_TO_WM2) -> float:
    return m_wm2 / met_wm2


def mifflin_bmr(sex: str, weight_kg: float, height_cm: float, age: float) -> float:
    base = 10 * weight_kg + 6.25 * height_cm - 5 * age
    return base + 5 if sex == "male" else base - 161


def du_bois_bsa(weight_kg: float, height_cm: float) -> float:
    return 0.007184 * (weight_kg ** 0.425) * (height_cm ** 0.725)


def estimate_metabolic_profile(
    sex: str,
    weight_kg: float,
    height_cm: float,
    age: float,
    rmr_factor: float = 1.1,
    activity_multiplier: float = 1.2,
    met_wm2: float = MET_TO_WM2,
) -> Dict[str, float]:
    bmr = mifflin_bmr(sex, weight_kg, height_cm, age)
    rmr = bmr * rmr_factor
    tee = rmr * activity_multiplier
    bsa = max(0.5, du_bois_bsa(weight_kg, height_cm))
    rmr_w = rmr * 0.0485
    tee_w = tee * 0.0485
    m_rmr = rmr_w / bsa
    m_active = tee_w / bsa
    return {
        "BMR_kcal_day": bmr,
        "RMR_kcal_day": rmr,
        "TEE_kcal_day": tee,
        "BSA_m2": bsa,
        "RMR_W": rmr_w,
        "TEE_W": tee_w,
        "M_rmr_Wm2": m_rmr,
        "M_active_Wm2": m_active,
        "met_active": wm2_to_met(m_active, met_wm2),
    }


def met_from_rmr(
    rmr_kcal_day: float,
    bsa_m2: float,
    activity_multiplier: float = 1.0,
    met_wm2: float = MET_TO_WM2,
) -> float:
    rmr_w = rmr_kcal_day * 0.0485
    m_wm2 = (rmr_w * activity_multiplier) / bsa_m2
    return wm2_to_met(m_wm2, met_wm2)
