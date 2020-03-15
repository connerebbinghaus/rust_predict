#[derive(Clone, Copy)]
pub struct EarthCenteredInertial {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Clone, Copy)]
pub struct ECIVector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl ECIVector {
    pub fn to_topo_horiz(
        self,
        ob: &super::latlong::LatitudeLongitude,
        now: chrono::DateTime<chrono::Utc>,
    ) -> super::th::TopographicHorizonVector {
        let theta = crate::time::gmst(now);
        let phi = ob.lat.to_radians();
        let (sin_theta, cos_theta) = theta.sin_cos();
        let (sin_phi, cos_phi) = phi.sin_cos();
        let s = sin_phi * cos_theta * self.x + sin_phi * sin_theta * self.y - cos_phi * self.z;
        let e = -sin_theta * self.x + cos_theta * self.y;
        let z = cos_phi * cos_theta * self.x + cos_phi * sin_theta * self.y + sin_phi * self.z;
        super::th::TopographicHorizonVector { s, e, z }
    }
}

impl EarthCenteredInertial {
    pub fn new(pos: (f64, f64, f64)) -> EarthCenteredInertial {
        EarthCenteredInertial {
            x: pos.0,
            y: pos.1,
            z: pos.2,
        }
    }
}

impl std::ops::Sub for EarthCenteredInertial {
    type Output = ECIVector;
    fn sub(self, other: EarthCenteredInertial) -> ECIVector {
        ECIVector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}
