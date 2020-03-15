/// Calculates the Julian date for a given point in time.
pub fn julian_date(date: chrono::DateTime<chrono::Utc>) -> f64 {
    use chrono::{Datelike, Timelike};

    let year = julian_date_of_year(date.year());
    let doy = julian_day_of_year(date.year(), date.month(), date.day());
    let f_day =
        julian_fraction_of_day(date.hour(), date.minute(), date.second(), date.nanosecond());
    year + doy as f64 + f_day
}

pub fn julian_date_of_year(year: i32) -> f64 {
    let year = year - 1;
    let a = year / 100;
    let b = 2 - a + (a / 4);
    (365.25 * year as f64).floor() + (30.6001f64 * 14.0).floor() + 1720994.5 + b as f64
}

const DAYS: [u16; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

pub fn julian_day_of_year(yr: i32, mo: u32, dy: u32) -> u16 {
    let mut day: u16 = 0;
    for i in 0..(mo - 1) as usize {
        day = day + DAYS[i];
    }

    let mut day = day + dy as u16;

    if ((yr % 4) == 0) && (((yr % 100) != 0) || ((yr % 400) == 0)) && (mo > 2) {
        day = day + 1;
    }

    day
}

pub fn julian_fraction_of_day(h: u32, m: u32, s: u32, ns: u32) -> f64 {
    ((h as f64 / 24.0)
        + (m as f64 / (24.0 * 60.0))
        + (s as f64 / (24.0 * 3600.0))
        + (ns as f64 / (24.0 * 3600.0 * 1.0e+9)))
    // - 0.5
}

/// Calculates the Greenwich Mean Sidereal Time (GMST) in radians at a point in time.
pub fn gmst(date: chrono::DateTime<chrono::Utc>) -> f64 {
    // Adapted from https://celestrak.com/c2olumns/v02n02/
    let jd = julian_date(date);
    let ut = (jd + 0.5).fract();
    let jd = jd - ut;
    let tu = (jd - 2451545.0) / 36525.0;
    let gmst = 24110.54841 + tu * (8640184.812866 + tu * (0.093104 - tu * 6.2E-6));
    let gmst = (gmst + 86400.0 * 1.00273790934 * ut).rem_euclid(86400.0);
    2.0 * std::f64::consts::PI * gmst / 86400.0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_julian_date() {
        use chrono::TimeZone;
        let date = chrono::Utc.ymd(1995, 10, 1).and_hms(9, 0, 0);
        let jd = julian_date(date);
        assert_approx_eq!(jd, 2449991.875);
    }

    #[test]
    fn test_gmst() {
        use chrono::TimeZone;
        let date = chrono::Utc.ymd(1995, 10, 1).and_hms(9, 0, 0);
        let gmst = gmst(date);
        assert_approx_eq!(gmst, 2.524218);
    }
}
