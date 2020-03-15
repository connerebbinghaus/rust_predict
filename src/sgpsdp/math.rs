use super::vals::*;

pub fn ac_tan(sinx: f64, cosx: f64) -> f64 {
    if cosx == 0.0 {
        if sinx > 0.0 {
            return PIO2;
        } else {
            return X3PIO2;
        }
    } else {
        if cosx > 0.0 {
            if sinx > 0.0 {
                return (sinx / cosx).atan();
            } else {
                return TWOPI + (sinx / cosx).atan();
            }
        } else {
            return PI + (sinx / cosx).atan();
        }
    }
}

/* Returns mod 2pi of argument */
pub fn fmod2p(x: f64) -> f64 {
    let i: i32;
    let mut ret_val: f64;

    ret_val = x;
    i = (ret_val / TWOPI) as i32;
    ret_val -= i as f64 * TWOPI;
    if ret_val < 0.0 {
        ret_val += TWOPI;
    }

    return ret_val;
}

/* The function ThetaG calculates the Greenwich Mean Sidereal Time */
/* for an epoch specified in the format used in the NORAD two-line */
/* element sets. It has now been adapted for dates beyond the year */
/* 1999, as described above. The function ThetaG_JD provides the   */
/* same calculation except that it is based on an input in the     */
/* form of a Julian Date. */
pub fn theta_g(epoch: f64, ds50: &mut f64) -> f64 {
    /* Reference:  The 1992 Astronomical Almanac, page B6. */

    /* Modification to support Y2K */
    /* Valid 1957 through 2056     */
    let mut year = (epoch * 1E-3).floor();
    let day = (epoch * 1E-3).fract() * 1E3;
    if year < 57.0 {
        year += 2000.0;
    } else {
        year += 1900.0;
    }
    /* End modification */

    let ut = day.fract();
    let day = day.floor();
    let jd = julian_date_of_year(year) + day;
    // tu = (jd - 2451545.0) / 36525.0;
    // gmst = 24110.54841 + tu * (8640184.812866 + tu * (0.093104 - tu * 6.2E-6));
    // gmst = (gmst + SECDAY * OMEGA_E * ut) % SECDAY;
    // thetag = TWOPI * gmst / SECDAY;
    *ds50 = jd - 2433281.5 + ut;

    fmod2p(6.3003880987 * *ds50 + 1.72944494)
}

/* The function Julian_Date_of_Year calculates the Julian Date  */
/* of Day 0.0 of {year}. This function is used to calculate the */
/* Julian Date of any date by using Julian_Date_of_Year, DOY,   */
/* and Fraction_of_Day. */
fn julian_date_of_year(mut year: f64) -> f64 {
    /* Astronomical Formulae for Calculators, Jean Meeus, */
    /* pages 23-25. Calculate Julian Date of 0.0 Jan year */

    year = year - 1.0;
    let mut i = (year / 100.0) as i64;
    let a = i;
    i = a / 4;
    let b = 2 - a + i;
    i = (365.25 * year) as i64;
    i += (30.6001 * 14.0) as i64;

    i as f64 + 1720994.5 + b as f64
}
