use crate::iapws97::constants;
// Region 5

const REGION_5_COEFFS_RES: [[f64; 3]; 7] = [
    [0.0, 0.0, 0.10658070028513e1],
    [1.0, 1.0, 0.15736404855259e-2],
    [1.0, 2.0, 0.90153761673944e-3],
    [1.0, 3.0, -0.50270077677648e-2],
    [2.0, 3.0, 0.22440037409485e-5],
    [2.0, 9.0, -0.41163275453471e-5],
    [3.0, 7.0, 0.37919454822955e-7],
];

const REGION_5_COEFFS_IDEAL: [[f64; 2]; 6] = [
    [0.0, -0.13179983674201e2],
    [1.0, 0.68540841634434e1],
    [-3.0, -0.24805148933466e-1],
    [-2.0, 0.36901534980333],
    [-1.0, -0.31161318213925e1],
    [2.0, -0.32961626538917],
];
// ================    Region 5 ===================

/// Returns the region-5 tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn tau_5(t: f64) -> f64 {
    1000.0 / t
}

/// Returns the region-2 pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn pi_5(p: f64) -> f64 {
    p / 1e6
}

///TODO: ideal functions look extremely similar to region 2
/// this can be extracted

/// Returns the region-2 ideal gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_5_ideal(t: f64, p: f64) -> f64 {
    let pi = pi_5(p);
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * tau.powf(ji);
    }
    pi.ln() + sum
}

/// Returns the region-5 ideal gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_5_ideal(_: f64, p: f64) -> f64 {
    1.0 / pi_5(p)
}

/// Returns the region-5 ideal gamma_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_5_ideal(t: f64, _: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * ji * tau.powf(ji - 1.0);
    }
    sum
}

/// Returns the region-5 ideal gamma_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_tau_5_ideal(t: f64, _: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    for coefficient in REGION_5_COEFFS_IDEAL {
        let ji = coefficient[0];
        let ni = coefficient[1];
        sum += ni * ji * (ji - 1.0) * tau.powf(ji - 2.0);
    }
    sum
}

/// Returns the region-5 residual gamma
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_5_res(t: f64, p: f64) -> f64 {
    let pi = pi_5(p);
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}

/// Returns the region-5 residual gamma_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_5_res(t: f64, p: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    let pi = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}

/// Returns the region-5 residual gamma_tau_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_tau_tau_5_res(t: f64, p: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    let pi = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}

/// Returns the region-5 residual gamma_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_5_res(t: f64, p: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    let pi = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}

/// Returns the region-5 residual gamma_pi_pi
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_pi_5_res(t: f64, p: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    let pi = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}

/// Returns the region-5 residual gamma_pi_tau
/// Temperature is assumed to be in K
/// Pressure is assumed to be in Pa
fn gamma_pi_tau_5_res(t: f64, p: f64) -> f64 {
    let mut sum: f64 = 0.0;
    let tau = tau_5(t);
    let pi = pi_5(p);
    for coefficient in REGION_5_COEFFS_RES {
        let ii = coefficient[0] as i32;
        let ji = coefficient[1] as i32;
        let ni = coefficient[2];
        sum += 0.0;
    }
    sum
}
