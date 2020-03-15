pub struct TopographicHorizonVector {
    pub s: f64,
    pub e: f64,
    pub z: f64,
}

impl TopographicHorizonVector {
    pub fn to_azel(self) -> super::azel::AzElVector {
        let range = (self.s * self.s + self.e * self.e + self.z * self.z).sqrt();
        let el = (self.z / range).asin();
        let mut az = (-self.e / self.s).atan();
        if self.s > 0.0 {
            az += std::f64::consts::PI;
        }
        if az < 0.0 {
            az += std::f64::consts::PI * 2.0;
        }
        super::azel::AzElVector { az, el, range }
    }
}
