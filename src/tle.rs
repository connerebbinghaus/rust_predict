use crate::sgpsdp::vals::*;

#[derive(Debug)]
pub struct TLE {
    name: String,
    sat_num: u32,
    class: char,
    desig: String,

    pub(crate) xno: f64,
    pub(crate) eo: f64,
    pub(crate) bstar: f64,
    pub(crate) xmo: f64,
    pub(crate) omegao: f64,
    pub(crate) xincl: f64,
    pub(crate) xnodeo: f64,
    pub(crate) omegao1: f64,
    pub(crate) xincl1: f64,
    pub(crate) xnodeo1: f64,
    pub(crate) epoch: f64,
    xndt2o: f64,
    xndd6o: f64,

    pub is_deep_space: bool,
    pub processed: bool,
}

impl TLE {
    pub fn parse(s: &str) -> TLE {
        tle::parse_tle(s).into()
    }

    pub fn epoch_jd(&self) -> f64 {
        let mut year = (self.epoch * 1E-3) as i32;
        let day = (self.epoch * 1E-3).fract() * 1E3;
        if year < 57 {
            year += 2000;
        } else {
            year += 1900;
        }
        /* End modification */

        crate::time::julian_date_of_year(year) + day
    }

    fn process_sgpsdp(&mut self) {
        let [ao, xnodp, mut dd1, dd2, delo, mut temp, a1, del1, r1]: [f64; 9];

        if self.processed {
            return;
        }

        /* Preprocess tle set */
        self.xnodeo *= DE2RA;
        self.omegao *= DE2RA;
        self.xmo *= DE2RA;
        self.xincl *= DE2RA;
        temp = TWOPI / XMNPDA / XMNPDA;

        self.xno = self.xno * temp * XMNPDA;
        self.xndt2o *= temp;
        self.xndd6o = self.xndd6o * temp / XMNPDA;
        self.bstar /= AE;

        /* Period > 225 minutes is deep space */
        dd1 = XKE / self.xno;
        dd2 = TOTHRD;
        a1 = dd1.powf(dd2);
        r1 = self.xincl.cos();
        dd1 = 1.0 - self.eo * self.eo;
        temp = CK2 * 1.5 * (r1 * r1 * 3.0 - 1.0) / dd1.powf(1.5);
        del1 = temp / (a1 * a1);
        ao = a1 * (1.0 - del1 * (TOTHRD * 0.5 + del1 * (del1 * 1.654320987654321 + 1.0)));
        delo = temp / (ao * ao);
        xnodp = self.xno / (delo + 1.0);

        /* Select a deep-space/near-earth ephemeris */
        if TWOPI / xnodp / XMNPDA >= 0.15625 {
            self.is_deep_space = true;
        } else {
            self.is_deep_space = false;
        }

        self.processed = true;
    }
}

impl From<tle::TLE> for TLE {
    fn from(tle: tle::TLE) -> Self {
        let mut r = TLE {
            name: tle.name,
            sat_num: tle.satellite_number,
            class: tle.classification,
            desig: tle.international_designator,
            epoch: tle.date.parse().unwrap(),
            xndt2o: tle.first_derivative_mean_motion,
            xndd6o: tle.second_derivative_mean_motion,
            bstar: tle.bstar,
            xincl: tle.inclination,
            xnodeo: tle.right_ascension,
            eo: tle.eccentricity,
            omegao: tle.argument_of_perigee,
            xmo: tle.mean_anomaly,
            xno: tle.mean_motion,
            omegao1: tle.argument_of_perigee,
            xincl1: tle.inclination,
            xnodeo1: tle.right_ascension,
            is_deep_space: false,
            processed: false,
        };
        r.process_sgpsdp();
        r
    }
}
