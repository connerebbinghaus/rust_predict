#[derive(Clone, Copy)]
pub struct LatitudeLongitude {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
}

impl LatitudeLongitude {
    pub fn to_eci(self, now: chrono::DateTime<chrono::Utc>) -> super::eci::EarthCenteredInertial {
        let pos = lla_to_eci(
            self.lat.to_radians(),
            self.long.to_radians(),
            self.alt / 1000.0,
            now,
        );
        super::eci::EarthCenteredInertial::new(pos)
    }
}

// Adapted from https://celestrak.com/columns/v02n02/
fn lla_to_eci(
    lat: f64,
    lon: f64,
    alt: f64,
    time: chrono::DateTime<chrono::Utc>,
) -> (f64, f64, f64) {
    const RE: f64 = 6378.135;
    let (sin_lat, cos_lat) = lat.sin_cos();
    let theta = (crate::time::gmst(time) + lon) % (2.0 * std::f64::consts::PI);
    let r = (RE + alt) * cos_lat;
    let (sin_theta, cos_theta) = theta.sin_cos();
    let x = r * cos_theta;
    let y = r * sin_theta;
    let z = (RE + alt) * sin_lat;
    (x, y, z)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::coords::EarthCenteredInertial;

    #[test]
    fn latlong_to_eci() {
        use chrono::TimeZone;
        let date = chrono::Utc.ymd(1995, 10, 1).and_hms(9, 0, 0);
        let lla = LatitudeLongitude {
            lat: 40.0,
            long: -75.0,
            alt: 0.0,
        };

        let EarthCenteredInertial { x, y, z } = lla.to_eci(date);
        assert_approx_eq!(x, 1700.938, 0.001);
        assert_approx_eq!(y, 4580.302, 0.001);
        assert_approx_eq!(z, 4099.786, 0.001);
    }
}
