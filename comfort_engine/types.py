from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, List


@dataclass(frozen=True)
class ComfortState:
    Ta: float
    Tr: float
    RH: float
    v: float
    patm: float = 101325.0


@dataclass(frozen=True)
class Occupant:
    m_wm2: float
    clo: float
    wme: float = 0.0
    posture: str = "seated"
    body_area_m2: Optional[float] = None


@dataclass(frozen=True)
class ModelConfig:
    standard: str = "ASHRAE_55"
    convection: str = "auto"
    hr_method: str = "auto"
    max_iter: int = 200
    tolerance: float = 0.00015


@dataclass(frozen=True)
class TraceStep:
    id: str
    equation: str
    inputs: Dict[str, float] = field(default_factory=dict)
    outputs: Dict[str, float] = field(default_factory=dict)
    units: Dict[str, str] = field(default_factory=dict)
    condition: Optional[str] = None


@dataclass
class PMVResult:
    PMV: float
    PPD: float
    Tcl: float
    h_c: float
    q_conv: float
    q_rad: float
    e_sk: float
    q_res: float
    L: float
    M: float
    H: float
    pmv_factor: float
    eq_slope: float
    dominance: str
    trace: List[TraceStep]
    warnings: List[str]
    branches: Dict[str, str]
    solver: Dict[str, float]
    components: Dict[str, float]
    components_w: Optional[Dict[str, float]] = None
    inputs: Dict[str, float] = field(default_factory=dict)
