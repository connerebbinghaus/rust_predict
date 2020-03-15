use crate::sgpsdp::Sat;
use crate::tle::TLE;

pub struct Satellite {
    tle: TLE,
    sgpsdp_data: Sat,
}

impl Satellite {
    pub fn from_tle(t: &str) -> Satellite {
        Satellite {
            tle: TLE::parse(t),
            sgpsdp_data: Sat::new(),
        }
    }

    pub fn predict(&mut self, date: chrono::DateTime<chrono::Utc>) {
        let dt = (crate::time::julian_date(date) - self.tle.epoch_jd()) * 24.0 * 60.0; // time since epoch in minutes
        self.sgpsdp_data.update(&mut self.tle, dt);
    }

    pub fn get_pos(&self) -> Option<crate::coords::EarthCenteredInertial> {
        self.sgpsdp_data.get_coords()
    }

    pub fn get_vel(&self) -> Option<crate::coords::eci::ECIVector> {
        self.sgpsdp_data.get_vel()
    }
}
