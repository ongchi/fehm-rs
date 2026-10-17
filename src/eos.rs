//! Liquid water density and viscosity: the rational polynomials of FEHM
//! (coefficient set 1 of `data.f`, evaluated as in `fluid_prop_control.f`).

/// Limits of the density fit of FEHM (coefficient set 1 of `data.f`)
pub const MIN_PRESSURE: f64 = 0.001;
pub const MAX_PRESSURE: f64 = 110.0;
pub const MIN_TEMPERATURE: f64 = 15.0;
pub const MAX_TEMPERATURE: f64 = 360.0;

/// Coefficients of one rational polynomial in pressure `p` (MPa) and
/// temperature `t` (°C), in the order of `crl`/`cvl` in `data.f`:
/// numerator `a0 + a1 p + a2 p² + a3 p³ + a4 t + a5 t² + a6 t³ + a7 p t + a8 p² t + a9 p t²`,
/// then the denominator with the same terms.
type RationalCoefficients = [f64; 20];

/// `crl(1:20, 1)` of `data.f`: liquid density (kg/m³)
const DENSITY: RationalCoefficients = [
    0.10000000e+01,
    0.17472599e-01,
    -0.20443098e-04,
    -0.17442012e-06,
    0.49564109e-02,
    -0.40757664e-04,
    0.50676664e-07,
    0.50330978e-04,
    0.33914814e-06,
    -0.18383009e-06,
    0.10009476e-02,
    0.16812589e-04,
    -0.24582622e-07,
    -0.17014984e-09,
    0.48841156e-05,
    -0.32967985e-07,
    0.28619380e-10,
    0.53249055e-07,
    0.30456698e-09,
    -0.12221899e-09,
];

/// `cvl(1:20, 1)` of `data.f`: liquid viscosity (Pa·s)
const VISCOSITY: RationalCoefficients = [
    0.17409149e-02,
    0.18894882e-04,
    -0.66439332e-07,
    -0.23122388e-09,
    -0.31534914e-05,
    0.11120716e-07,
    -0.48576020e-10,
    0.28006861e-07,
    0.23225035e-09,
    0.47180171e-10,
    0.10000000e+01,
    0.10523153e-01,
    -0.22658391e-05,
    -0.31796607e-06,
    0.29869141e-01,
    0.21844248e-03,
    -0.87658855e-06,
    0.41690362e-03,
    -0.25147022e-05,
    0.22144660e-05,
];

/// A property and its derivative with respect to pressure (per MPa)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Property {
    pub value: f64,
    pub pressure_derivative: f64,
}

/// Liquid water density (kg/m³) at pressure `pressure` (MPa) and temperature `temperature` (°C)
pub fn density(pressure: f64, temperature: f64) -> Property {
    evaluate(&DENSITY, pressure, temperature)
}

/// Liquid water viscosity (Pa·s) at pressure `pressure` (MPa) and temperature `temperature` (°C)
pub fn viscosity(pressure: f64, temperature: f64) -> Property {
    evaluate(&VISCOSITY, pressure, temperature)
}

/// Value and pressure derivative of a rational polynomial, with the same
/// order of operations as `fluid_prop_control.f`
fn evaluate(c: &RationalCoefficients, p: f64, t: f64) -> Property {
    let p2 = p * p;
    let p3 = p2 * p;
    let t2 = t * t;
    let t3 = t2 * t;
    let tp = p * t;
    let t2p = t2 * p;
    let tp2 = t * p2;

    let numerator = (c[0] + c[1] * p + c[2] * p2 + c[3] * p3)
        + (c[4] * t + c[5] * t2 + c[6] * t3)
        + (c[7] * tp + c[9] * t2p + c[8] * tp2);
    let denominator = (c[10] + c[11] * p + c[12] * p2 + c[13] * p3)
        + (c[14] * t + c[15] * t2 + c[16] * t3)
        + (c[17] * tp + c[19] * t2p + c[18] * tp2);

    let numerator_derivative =
        c[1] + 2.0 * c[2] * p + 3.0 * c[3] * p2 + c[7] * t + 2.0 * c[8] * tp + c[9] * t2;
    let denominator_derivative =
        c[11] + 2.0 * c[12] * p + 3.0 * c[13] * p2 + c[17] * t + 2.0 * c[18] * tp + c[19] * t2;

    Property {
        value: numerator / denominator,
        pressure_derivative: (denominator * numerator_derivative
            - numerator * denominator_derivative)
            / (denominator * denominator),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_derivative_matches_finite_difference(property: fn(f64, f64) -> Property) {
        let step = 1e-4;
        for pressure in [0.1, 1.0, 10.0, 50.0, 100.0] {
            for temperature in [20.0, 50.0, 150.0, 300.0] {
                let finite_difference = (property(pressure + step, temperature).value
                    - property(pressure - step, temperature).value)
                    / (2.0 * step);
                let derivative = property(pressure, temperature).pressure_derivative;
                assert!(
                    (derivative - finite_difference).abs() <= 1e-6 * finite_difference.abs(),
                    "P = {pressure}, T = {temperature}: {derivative} != {finite_difference}"
                );
            }
        }
    }

    #[test]
    fn density_derivative_matches_finite_difference() {
        assert_derivative_matches_finite_difference(density);
    }

    #[test]
    fn viscosity_derivative_matches_finite_difference() {
        assert_derivative_matches_finite_difference(viscosity);
    }

    #[test]
    fn density_and_viscosity_are_near_tabulated_water_values() {
        // NIST: 10 MPa, 50 °C: 992.6 kg/m³, 5.49e-4 Pa·s; the fit has an error of up to 0.34 %
        let rho = density(10.0, 50.0).value;
        let mu = viscosity(10.0, 50.0).value;
        assert!((rho - 992.6).abs() < 0.01 * 992.6, "density {rho}");
        assert!((mu - 5.49e-4).abs() < 0.01 * 5.49e-4, "viscosity {mu}");
    }
}
