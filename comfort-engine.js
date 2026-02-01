export const MET_TO_WM2 = 58.15;
export const CLO_TO_ICL = 0.155;
export const SIGMA = 5.67e-8;
export const EPSILON = 0.95;
export const F_EFF = 0.72;

export const metToWm2 = (met) => met * MET_TO_WM2;
export const wm2ToMet = (mWm2) => mWm2 / MET_TO_WM2;

export const mifflinBMR = ({ sex, weightKg, heightCm, age }) => {
  const base = 10 * weightKg + 6.25 * heightCm - 5 * age;
  return sex === 'male' ? base + 5 : base - 161;
};

export const duBoisBSA = ({ weightKg, heightCm }) => 0.007184 * Math.pow(weightKg, 0.425) * Math.pow(heightCm, 0.725);

export const estimateMetabolicProfile = ({
  sex,
  weightKg,
  heightCm,
  age,
  rmrFactor = 1.1,
  activityMultiplier = 1.2
}) => {
  const BMR = mifflinBMR({ sex, weightKg, heightCm, age });
  const RMR = BMR * rmrFactor;
  const TEE = RMR * activityMultiplier;
  const BSA = Math.max(0.5, duBoisBSA({ weightKg, heightCm }));
  const RMR_W = RMR * 0.0485;
  const TEE_W = TEE * 0.0485;
  const M_rmr = RMR_W / BSA;
  const M_active = TEE_W / BSA;
  return {
    BMR,
    RMR,
    TEE,
    BSA,
    RMR_W,
    TEE_W,
    M_rmr,
    M_active,
    met_active: M_active / MET_TO_WM2
  };
};

export const satVaporPressure = (taC) => 610.94 * Math.exp((17.625 * taC) / (taC + 243.04));

export const dewPointFromRh = (taC, rhPct) => {
  if (!Number.isFinite(rhPct) || rhPct <= 0) return NaN;
  const gamma = Math.log(rhPct / 100) + (17.625 * taC) / (243.04 + taC);
  return (243.04 * gamma) / (17.625 - gamma);
};

export const humidityRatio = (p_a, p_atm = 101325) => 0.62198 * p_a / (p_atm - p_a);

export const computePMV = (env, occ) => {
  const Ta = env.Ta;
  const Tr = env.Tr;
  const RH = env.RH;
  const v = Number.isFinite(env.v) ? env.v : 0;
  const patm = env.patm ?? 101325;

  const m_wm2 = occ.m_wm2;
  const clo = occ.clo;
  const wme = occ.wme ?? 0;

  if (!Number.isFinite(m_wm2)) {
    throw new Error('computePMV expects occ.m_wm2 in W/m².');
  }

  const M = m_wm2;
  const H = M - wme;

  const Icl = clo * CLO_TO_ICL;
  const f_cl = Icl <= 0.078 ? 1 + 1.29 * Icl : 1.05 + 0.645 * Icl;

  const p_ws = satVaporPressure(Ta);
  const p_a = (RH / 100) * p_ws;

  const h_forced = 12.1 * Math.sqrt(Math.max(v, 0));
  const taa = Ta + 273;
  const tra = Tr + 273;
  const tcla = taa + (35.5 - Ta) / (3.5 * Icl + 0.1);

  const p1 = Icl * f_cl;
  const p2 = p1 * 3.96;
  const p3 = p1 * 100;
  const p4 = p1 * taa;
  const p5 = 308.7 - 0.028 * H + p2 * Math.pow(tra / 100, 4);

  let xn = tcla / 100;
  let xf = xn;
  let h_c = h_forced;
  let iterations = 0;
  let residual = Infinity;
  for (let i = 0; i < 200; i++) {
    xf = (xf + xn) / 2;
    const h_natural_iter = 2.38 * Math.pow(Math.abs(100 * xf - taa), 0.25);
    h_c = Math.max(h_forced, h_natural_iter);
    xn = (p5 + p4 * h_c - p2 * Math.pow(xf, 4)) / (100 + p3 * h_c);
    residual = Math.abs(xn - xf);
    iterations = i + 1;
    if (residual <= 0.00015) {
      break;
    }
  }

  const Tcl = 100 * xn - 273;
  const h_natural = 2.38 * Math.pow(Math.abs(Tcl - Ta), 0.25);
  h_c = Math.max(h_forced, h_natural);

  const Tcl_K = Tcl + 273.15;
  const Tr_K = Tr + 273.15;
  const Q_rad = 3.96e-8 * f_cl * (Math.pow(Tcl_K, 4) - Math.pow(Tr_K, 4));
  const Q_conv = f_cl * h_c * (Tcl - Ta);

  const E_diff = 3.05e-3 * (5733 - 6.99 * H - p_a);
  const E_sw = Math.max(0, 0.42 * (H - MET_TO_WM2));
  const E_res = 1.7e-5 * M * (5867 - p_a);
  const C_res = 0.0014 * M * (34 - Ta);

  const E_sk = E_diff + E_sw;
  const Q_res = E_res + C_res;
  const L = Q_conv + Q_rad + E_sk + Q_res;

  const pmvFactor = 0.303 * Math.exp(-0.036 * M) + 0.028;
  const PMV = pmvFactor * (H - L);
  const PPD = 100 - 95 * Math.exp(-0.03353 * Math.pow(PMV, 4) - 0.2179 * Math.pow(PMV, 2));

  const deltaTr = Tcl - Tr;
  let hr = NaN;
  let hrMethod = 'direct';
  if (Number.isFinite(deltaTr) && Math.abs(deltaTr) > 0.05) {
    hr = Q_rad / (f_cl * deltaTr);
  } else {
    const tMean = (Tcl + Tr) / 2 + 273.15;
    hr = 4 * EPSILON * SIGMA * F_EFF * Math.pow(tMean, 3);
    hrMethod = 'linearized';
  }

  const eqSlope = Number.isFinite(hr) && Number.isFinite(h_c) ? -hr / h_c : NaN;
  let dominance = '—';
  if (Number.isFinite(hr) && Number.isFinite(h_c) && h_c > 0) {
    const ratio = hr / h_c;
    if (ratio > 1.2) dominance = 'Radiation-dominant';
    else if (ratio < 0.8) dominance = 'Convection-dominant';
    else dominance = 'Mixed';
  }

  return {
    M,
    H,
    Icl,
    f_cl,
    p_ws,
    p_a,
    Tcl,
    h_c,
    h_forced,
    h_natural,
    Q_conv,
    Q_rad,
    E_diff,
    E_sw,
    E_res,
    C_res,
    E_sk,
    Q_res,
    L,
    PMV,
    PPD,
    pmvFactor,
    hr,
    eqSlope,
    dominance,
    solver: {
      iterations,
      residual,
      converged: residual <= 0.00015
    },
    branches: {
      convectionRegime: h_forced >= h_natural ? 'forced' : 'natural',
      hrMethod
    },
    patm
  };
};
