/*
 *  Unit SGP4SDP4
 *           Author:  Dr TS Kelso
 * Original Version:  1991 Oct 30
 * Current Revision:  1992 Sep 03
 *          Version:  1.50
 *        Copyright:  1991-1992, All Rights Reserved
 *
 *   Ported to C by:  Neoklis Kyriazis  April 10  2001
 *   Reentrancy mods by Alexandru Csete OZ9AEC
 *   Rust port by Conner Ebbinghaus
 */

mod math;
#[cfg(test)]
mod tests;
pub(crate) mod vals;

use math::*;
use vals::*;

use crate::tle::TLE;

pub struct Sat {
    data: SDPSGPData,
}

impl Sat {
    pub fn new() -> Sat {
        Sat {
            data: SDPSGPData::None,
        }
    }

    pub fn update(&mut self, tle: &mut crate::tle::TLE, dt: f64) {
        if tle.is_deep_space {
            sdp4(&mut self.data, tle, dt);
        } else {
            sgp4(&mut self.data, tle, dt);
        }
    }

    pub fn get_coords(&self) -> Option<crate::coords::EarthCenteredInertial> {
        let pos = match &self.data {
            SDPSGPData::SGP(sgp) => sgp.pos,
            SDPSGPData::SDP { sdp, .. } => sdp.pos,
            SDPSGPData::None => return None,
        };

        Some(crate::coords::EarthCenteredInertial::new((
            pos.0 * XKMPER,
            pos.1 * XKMPER,
            pos.2 * XKMPER,
        )))
    }

    pub fn get_vel(&self) -> Option<crate::coords::eci::ECIVector> {
        let vel = match &self.data {
            SDPSGPData::SGP(sgp) => sgp.vel,
            SDPSGPData::SDP { sdp, .. } => sdp.vel,
            SDPSGPData::None => return None,
        };

        Some(crate::coords::eci::ECIVector {
            x: vel.0 * (XKMPER * XMNPDA / SECDAY),
            y: vel.1 * (XKMPER * XMNPDA / SECDAY),
            z: vel.2 * (XKMPER * XMNPDA / SECDAY),
        })
    }
}

enum SDPSGPData {
    None,
    SGP(SGPData),
    SDP { sdp: SDPData, dps: DeepData },
}

impl SDPSGPData {
    fn is_sgp(&self) -> bool {
        match self {
            SDPSGPData::SGP(_) => true,
            _ => false,
        }
    }

    fn assert_sgp(&mut self) -> &mut SGPData {
        match self {
            SDPSGPData::SGP(ref mut d) => d,
            _ => unreachable!(),
        }
    }

    fn is_sdp(&self) -> bool {
        match self {
            SDPSGPData::SDP { .. } => true,
            _ => false,
        }
    }

    fn assert_sdp(&mut self) -> (&mut SDPData, &mut DeepData) {
        match self {
            SDPSGPData::SDP { sdp, dps } => (sdp, dps),
            _ => unreachable!(),
        }
    }
}

#[derive(Default)]
struct SGPData {
    simple: bool,
    cosio: f64,
    x3thm1: f64,
    xnodp: f64,
    aodp: f64,
    eta: f64,
    c1: f64,
    c4: f64,
    c5: f64,
    sinio: f64,
    x1mth2: f64,
    xmdot: f64,
    omgdot: f64,
    xnodot: f64,
    omgcof: f64,
    xmcof: f64,
    xnodcf: f64,
    t2cof: f64,
    xlcof: f64,
    aycof: f64,
    delmo: f64,
    sinmo: f64,
    x7thm1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    pos: (f64, f64, f64),
    vel: (f64, f64, f64),
    phase: f64,
}

/* SGP4 */
/* This function is used to calculate the position and velocity */
/* of near-earth (period < 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s.*/
fn sgp4(data: &mut SDPSGPData, tle: &mut TLE, tsince: f64) {
    /* Initialization */
    let sgp = if !data.is_sgp() {
        let mut sgp = SGPData::default();

        /* Recover original mean motion (xnodp) and   */
        /* semimajor axis (aodp) from input elements. */
        let a1 = (XKE / tle.xno).powf(TOTHRD);
        sgp.cosio = tle.xincl.cos();
        let theta2 = sgp.cosio * sgp.cosio;
        sgp.x3thm1 = 3.0 * theta2 - 1.0;
        let eosq = tle.eo * tle.eo;
        let betao2 = 1.0 - eosq;
        let betao = betao2.sqrt();
        let del1 = 1.5 * CK2 * sgp.x3thm1 / (a1 * a1 * betao * betao2);
        let ao = a1 * (1.0 - del1 * (0.5 * TOTHRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
        let delo = 1.5 * CK2 * sgp.x3thm1 / (ao * ao * betao * betao2);
        sgp.xnodp = tle.xno / (1.0 + delo);
        sgp.aodp = ao / (1.0 - delo);

        /* For perigee less than 220 kilometers, the "simple" flag is set */
        /* and the equations are truncated to linear variation in sqrt a  */
        /* and quadratic variation in mean anomaly.  Also, the c3 term,   */
        /* the delta omega term, and the delta m term are dropped.        */
        if (sgp.aodp * (1.0 - tle.eo) / AE) < (220.0 / XKMPER + AE) {
            sgp.simple = true;
        } else {
            sgp.simple = false;
        }

        /* For perigee below 156 km, the       */
        /* values of s and qoms2t are altered. */
        let mut s4 = __S__;
        let mut qoms24 = QOMS2T;
        let perige = (sgp.aodp * (1.0 - tle.eo) - AE) * XKMPER;
        if perige < 156.0 {
            if perige <= 98.0 {
                s4 = 20.0;
            } else {
                s4 = perige - 78.0;
            }
            qoms24 = ((120.0 - s4) * AE / XKMPER).powi(4);
            s4 = s4 / XKMPER + AE;
        };

        let pinvsq = 1.0 / (sgp.aodp * sgp.aodp * betao2 * betao2);
        let tsi = 1.0 / (sgp.aodp - s4);
        sgp.eta = sgp.aodp * tle.eo * tsi;
        let etasq = sgp.eta * sgp.eta;
        let eeta = tle.eo * sgp.eta;
        let psisq = (1.0 - etasq).abs();
        let coef = qoms24 * tsi.powi(4);
        let coef1 = coef / psisq.powf(3.5);
        let c2 = coef1
            * sgp.xnodp
            * (sgp.aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                + 0.75 * CK2 * tsi / psisq * sgp.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        sgp.c1 = c2 * tle.bstar;
        sgp.sinio = (tle.xincl).sin();
        let a3ovk2 = -XJ3 / CK2 * AE.powi(3);
        let c3 = coef * tsi * a3ovk2 * sgp.xnodp * AE * sgp.sinio / tle.eo;
        sgp.x1mth2 = 1.0 - theta2;
        sgp.c4 = 2.0
            * sgp.xnodp
            * coef1
            * sgp.aodp
            * betao2
            * (sgp.eta * (2.0 + 0.5 * etasq) + tle.eo * (0.5 + 2.0 * etasq)
                - 2.0 * CK2 * tsi / (sgp.aodp * psisq)
                    * (-3.0 * sgp.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                        + 0.75
                            * sgp.x1mth2
                            * (2.0 * etasq - eeta * (1.0 + etasq))
                            * (2.0 * tle.omegao).cos()));
        sgp.c5 = 2.0 * coef1 * sgp.aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        let theta4 = theta2 * theta2;
        let temp1 = 3.0 * CK2 * pinvsq * sgp.xnodp;
        let temp2 = temp1 * CK2 * pinvsq;
        let temp3 = 1.25 * CK4 * pinvsq * pinvsq * sgp.xnodp;
        sgp.xmdot = sgp.xnodp
            + 0.5 * temp1 * betao * sgp.x3thm1
            + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
        let x1m5th = 1.0 - 5.0 * theta2;
        sgp.omgdot = -0.5 * temp1 * x1m5th
            + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4)
            + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
        let xhdot1 = -temp1 * sgp.cosio;
        sgp.xnodot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2))
                * sgp.cosio;
        sgp.omgcof = tle.bstar * c3 * tle.omegao.cos();
        sgp.xmcof = -TOTHRD * coef * tle.bstar * AE / eeta;
        sgp.xnodcf = 3.5 * betao2 * xhdot1 * sgp.c1;
        sgp.t2cof = 1.5 * sgp.c1;
        sgp.xlcof = 0.125 * a3ovk2 * sgp.sinio * (3.0 + 5.0 * sgp.cosio) / (1.0 + sgp.cosio);
        sgp.aycof = 0.25 * a3ovk2 * sgp.sinio;
        sgp.delmo = (1.0 + sgp.eta * tle.xmo.cos()).powi(3);
        sgp.sinmo = tle.xmo.sin();
        sgp.x7thm1 = 7.0 * theta2 - 1.0;
        if sgp.simple {
            let c1sq = sgp.c1 * sgp.c1;
            sgp.d2 = 4.0 * sgp.aodp * tsi * c1sq;
            let temp = sgp.d2 * tsi * sgp.c1 / 3.0;
            sgp.d3 = (17.0 * sgp.aodp + s4) * temp;
            sgp.d4 = 0.5 * temp * sgp.aodp * tsi * (221.0 * sgp.aodp + 31.0 * s4) * sgp.c1;
            sgp.t3cof = sgp.d2 + 2.0 * c1sq;
            sgp.t4cof = 0.25 * (3.0 * sgp.d3 + sgp.c1 * (12.0 * sgp.d2 + 10.0 * c1sq));
            sgp.t5cof = 0.2
                * (3.0 * sgp.d4
                    + 12.0 * sgp.c1 * sgp.d3
                    + 6.0 * sgp.d2 * sgp.d2
                    + 15.0 * c1sq * (2.0 * sgp.d2 + c1sq));
        };

        *data = SDPSGPData::SGP(sgp);
        data.assert_sgp()
    } else {
        data.assert_sgp()
    };

    /* Update for secular gravity and atmospheric drag. */
    let xmdf = tle.xmo + sgp.xmdot * tsince;
    let omgadf = tle.omegao + sgp.omgdot * tsince;
    let xnoddf = tle.xnodeo + sgp.xnodot * tsince;
    let mut omega = omgadf;
    let mut xmp = xmdf;
    let tsq = tsince * tsince;
    let xnode = xnoddf + sgp.xnodcf * tsq;
    let mut tempa = 1.0 - sgp.c1 * tsince;
    let mut tempe = tle.bstar * sgp.c4 * tsince;
    let mut templ = sgp.t2cof * tsq;
    if sgp.simple {
        let delomg = sgp.omgcof * tsince;
        let delm = sgp.xmcof * ((1.0 + sgp.eta * xmdf.cos()).powi(3) - sgp.delmo);
        let temp = delomg + delm;
        xmp = xmdf + temp;
        omega = omgadf - temp;
        let tcube = tsq * tsince;
        let tfour = tsince * tcube;
        tempa = tempa - sgp.d2 * tsq - sgp.d3 * tcube - sgp.d4 * tfour;
        tempe = tempe + tle.bstar * sgp.c5 * (xmp.sin() - sgp.sinmo);
        templ = templ + sgp.t3cof * tcube + tfour * (sgp.t4cof + tsince * sgp.t5cof);
    };

    let a = sgp.aodp * tempa.powi(2);
    let e = tle.eo - tempe;
    let xl = xmp + omega + xnode + sgp.xnodp * templ;
    let beta = (1.0 - e * e).sqrt();
    let xn = XKE / a.powf(1.5);

    /* Long period periodics */
    let axn = e * omega.cos();
    let temp = 1.0 / (a * beta * beta);
    let xll = temp * sgp.xlcof * axn;
    let aynl = temp * sgp.aycof;
    let xlt = xl + xll;
    let ayn = e * omega.sin() + aynl;

    /* Solve Kepler's' Equation */
    let capu = fmod2p(xlt - xnode);
    let mut temp2 = capu;

    let mut temp3 = 0.0;
    let mut temp4 = 0.0;
    let mut temp5 = 0.0;
    let mut temp6 = 0.0;
    let mut sinepw= 0.0;
    let mut cosepw= 0.0;

    for _i in 0..=10 {
        let (sinepw_t, cosepw_t) = temp2.sin_cos();
        sinepw = sinepw_t;
        cosepw = cosepw_t;
        temp3 = axn * sinepw;
        temp4 = ayn * cosepw;
        temp5 = axn * cosepw;
        temp6 = ayn * sinepw;
        let epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
        if (epw - temp2).abs() <= E6A {
            break;
        }
        temp2 = epw;
    }

    /* Short period preliminary quantities */
    let ecose = temp5 + temp6;
    let esine = temp3 - temp4;
    let elsq = axn * axn + ayn * ayn;
    let temp = 1.0 - elsq;
    let pl = a * temp;
    let r = a * (1.0 - ecose);
    let temp1 = 1.0 / r;
    let rdot = XKE * a.sqrt() * esine * temp1;
    let rfdot = XKE * pl.sqrt() * temp1;
    temp2 = a * temp1;
    let betal = temp.sqrt();
    temp3 = 1.0 / (1.0 + betal);
    let cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
    let sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
    let u = ac_tan(sinu, cosu);
    let sin2u = 2.0 * sinu * cosu;
    let cos2u = 2.0 * cosu * cosu - 1.0;
    let temp = 1.0 / pl;
    let temp1 = CK2 * temp;
    temp2 = temp1 * temp;

    /* Update for short periodics */
    let rk = r * (1.0 - 1.5 * temp2 * betal * sgp.x3thm1) + 0.5 * temp1 * sgp.x1mth2 * cos2u;
    let uk = u - 0.25 * temp2 * sgp.x7thm1 * sin2u;
    let xnodek = xnode + 1.5 * temp2 * sgp.cosio * sin2u;
    let xinck = tle.xincl + 1.5 * temp2 * sgp.cosio * sgp.sinio * cos2u;
    let rdotk = rdot - xn * temp1 * sgp.x1mth2 * sin2u;
    let rfdotk = rfdot + xn * temp1 * (sgp.x1mth2 * cos2u + 1.5 * sgp.x3thm1);

    /* Orientation vectors */
    let (sinuk, cosuk) = uk.sin_cos();
    let (sinik, cosik) = xinck.sin_cos();
    let (sinnok, cosnok) = xnodek.sin_cos();
    let xmx = -sinnok * cosik;
    let xmy = cosnok * cosik;
    let ux = xmx * sinuk + cosnok * cosuk;
    let uy = xmy * sinuk + sinnok * cosuk;
    let uz = sinik * sinuk;
    let vx = xmx * cosuk - cosnok * sinuk;
    let vy = xmy * cosuk - sinnok * sinuk;
    let vz = sinik * cosuk;

    /* Position and velocity */
    sgp.pos.0 = rk * ux;
    sgp.pos.1 = rk * uy;
    sgp.pos.2 = rk * uz;
    sgp.vel.0 = rdotk * ux + rfdotk * vx;
    sgp.vel.1 = rdotk * uy + rfdotk * vy;
    sgp.vel.2 = rdotk * uz + rfdotk * vz;

    sgp.phase = xlt - xnode - omgadf + TWOPI;
    if sgp.phase < 0.0 {
        sgp.phase += TWOPI;
    }
    sgp.phase = fmod2p(sgp.phase);

    tle.omegao1 = omega;
    tle.xincl1 = xinck;
    tle.xnodeo1 = xnodek;
}

#[derive(Default)]
struct SDPData {
    cosio: f64,
    theta2: f64,
    x3thm1: f64,
    eosq: f64,
    betao2: f64,
    betao: f64,
    xnodp: f64,
    aodp: f64,
    sing: f64,
    cosg: f64,
    c1: f64,
    c4: f64,
    sinio: f64,
    x1mth2: f64,
    xmdot: f64,
    omgdot: f64,
    xnodot: f64,
    xnodcf: f64,
    t2cof: f64,
    xlcof: f64,
    aycof: f64,
    x7thm1: f64,
    omgadf: f64,
    xnode: f64,
    xn: f64,
    xll: f64,
    t: f64,
    em: f64,
    xinc: f64,
    ds50: f64,
    pos: (f64, f64, f64),
    vel: (f64, f64, f64),
    phase: f64,
}

/* SDP4 */
/* This function is used to calculate the position and velocity */
/* of deep-space (period > 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s. */
fn sdp4(data: &mut SDPSGPData, tle: &mut TLE, tsince: f64) {


    /* Initialization */
    let (mut sdp, dps) = if !data.is_sdp() {
        let mut sdp = SDPData::default();
        let mut dps = DeepData::default();

        /* Recover original mean motion (xnodp) and   */
        /* semimajor axis (aodp) from input elements. */
        let a1 = (XKE / tle.xno).powf(TOTHRD);
        sdp.cosio = tle.xincl.cos();
        sdp.theta2 = sdp.cosio * sdp.cosio;
        sdp.x3thm1 = 3.0 * sdp.theta2 - 1.0;
        sdp.eosq = tle.eo * tle.eo;
        sdp.betao2 = 1.0 - sdp.eosq;
        sdp.betao = sdp.betao2.sqrt();
        let del1 = 1.5 * CK2 * sdp.x3thm1 / (a1 * a1 * sdp.betao * sdp.betao2);
        let ao = a1 * (1.0 - del1 * (0.5 * TOTHRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
        let delo = 1.5 * CK2 * sdp.x3thm1 / (ao * ao * sdp.betao * sdp.betao2);
        sdp.xnodp = tle.xno / (1.0 + delo);
        sdp.aodp = ao / (1.0 - delo);

        /* For perigee below 156 km, the values */
        /* of s and qoms2t are altered.         */
        let mut s4 = __S__;
        let mut qoms24 = QOMS2T;
        let perige = (sdp.aodp * (1.0 - tle.eo) - AE) * XKMPER;
        if perige < 156.0 {
            if perige <= 98.0 {
                s4 = 20.0;
            } else {
                s4 = perige - 78.0;
            }
            qoms24 = ((120.0 - s4) * AE / XKMPER).powi(4);
            s4 = s4 / XKMPER + AE;
        }
        let pinvsq = 1.0 / (sdp.aodp * sdp.aodp * sdp.betao2 * sdp.betao2);
        let (sing_t, cosg_t) = tle.omegao.sin_cos();
        sdp.sing = sing_t;
        sdp.cosg = cosg_t;
        let tsi = 1.0 / (sdp.aodp - s4);
        let eta = sdp.aodp * tle.eo * tsi;
        let etasq = eta * eta;
        let eeta = tle.eo * eta;
        let psisq = (1.0 - etasq).abs();
        let coef = qoms24 * tsi.powi(4);
        let coef1 = coef / psisq.powf(3.5);
        let c2 = coef1
            * sdp.xnodp
            * (sdp.aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                + 0.75 * CK2 * tsi / psisq * sdp.x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
        sdp.c1 = tle.bstar * c2;
        sdp.sinio = tle.xincl.sin();
        let a3ovk2 = -XJ3 / CK2 * AE.powi(3);
        sdp.x1mth2 = 1.0 - sdp.theta2;
        sdp.c4 = 2.0
            * sdp.xnodp
            * coef1
            * sdp.aodp
            * sdp.betao2
            * (eta * (2.0 + 0.5 * etasq) + tle.eo * (0.5 + 2.0 * etasq)
                - 2.0 * CK2 * tsi / (sdp.aodp * psisq)
                    * (-3.0 * sdp.x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                        + 0.75
                            * sdp.x1mth2
                            * (2.0 * etasq - eeta * (1.0 + etasq))
                            * (2.0 * tle.omegao).cos()));
        let theta4 = sdp.theta2 * sdp.theta2;
        let temp1 = 3.0 * CK2 * pinvsq * sdp.xnodp;
        let temp2 = temp1 * CK2 * pinvsq;
        let temp3 = 1.25 * CK4 * pinvsq * pinvsq * sdp.xnodp;
        sdp.xmdot = sdp.xnodp
            + 0.5 * temp1 * sdp.betao * sdp.x3thm1
            + 0.0625 * temp2 * sdp.betao * (13.0 - 78.0 * sdp.theta2 + 137.0 * theta4);
        let x1m5th = 1.0 - 5.0 * sdp.theta2;
        sdp.omgdot = -0.5 * temp1 * x1m5th
            + 0.0625 * temp2 * (7.0 - 114.0 * sdp.theta2 + 395.0 * theta4)
            + temp3 * (3.0 - 36.0 * sdp.theta2 + 49.0 * theta4);
        let xhdot1 = -temp1 * sdp.cosio;
        sdp.xnodot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * sdp.theta2) + 2.0 * temp3 * (3.0 - 7.0 * sdp.theta2))
                * sdp.cosio;
        sdp.xnodcf = 3.5 * sdp.betao2 * xhdot1 * sdp.c1;
        sdp.t2cof = 1.5 * sdp.c1;
        sdp.xlcof = 0.125 * a3ovk2 * sdp.sinio * (3.0 + 5.0 * sdp.cosio) / (1.0 + sdp.cosio);
        sdp.aycof = 0.25 * a3ovk2 * sdp.sinio;
        sdp.x7thm1 = 7.0 * sdp.theta2 - 1.0;

        /* initialize Deep() */
        deep(DPEntry::Init, &mut sdp, &mut dps, tle);
        *data = SDPSGPData::SDP { sdp, dps };
        data.assert_sdp()
    } else {
        data.assert_sdp()
    };

    /* Update for secular gravity and atmospheric drag */
    let xmdf = tle.xmo + sdp.xmdot * tsince;
    sdp.omgadf = tle.omegao + sdp.omgdot * tsince;
    let xnoddf = tle.xnodeo + sdp.xnodot * tsince;
    let tsq = tsince * tsince;
    sdp.xnode = xnoddf + sdp.xnodcf * tsq;
    let tempa = 1.0 - sdp.c1 * tsince;
    let tempe = tle.bstar * sdp.c4 * tsince;
    let templ = sdp.t2cof * tsq;
    sdp.xn = sdp.xnodp;

    /* Update for deep-space secular effects */
    sdp.xll = xmdf;
    sdp.t = tsince;

    deep(DPEntry::Sec, sdp, dps, tle);

    let xmdf = sdp.xll;
    let a = (XKE / sdp.xn).powf(TOTHRD) * tempa * tempa;
    sdp.em = sdp.em - tempe;
    let xmam = xmdf + sdp.xnodp * templ;

    /* Update for deep-space periodic effects */
    sdp.xll = xmam;

    deep(DPEntry::Per, sdp, dps, tle);

    let xmam = sdp.xll;
    let xl = xmam + sdp.omgadf + sdp.xnode;
    let beta = (1.0 - sdp.em * sdp.em).sqrt();
    sdp.xn = XKE / a.powf(1.5);

    /* Long period periodics */
    let axn = sdp.em * sdp.omgadf.cos();
    let temp = 1.0 / (a * beta * beta);
    let xll = temp * sdp.xlcof * axn;
    let aynl = temp * sdp.aycof;
    let xlt = xl + xll;
    let ayn = sdp.em * sdp.omgadf.sin() + aynl;

    /* Solve Kepler's Equation */
    let capu = fmod2p(xlt - sdp.xnode);
    let mut temp2 = capu;

    let mut temp3 = 0.0;
    let mut temp4 = 0.0;
    let mut temp5 = 0.0;
    let mut temp6 = 0.0;
    let mut sinepw= 0.0;
    let mut cosepw= 0.0;

    for _i in 0..=10 {
        let (sinepw_t, cosepw_t) = temp2.sin_cos();
        sinepw = sinepw_t;
        cosepw = cosepw_t;
        temp3 = axn * sinepw;
        temp4 = ayn * cosepw;
        temp5 = axn * cosepw;
        temp6 = ayn * sinepw;
        let epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
        if (epw - temp2).abs() <= E6A {
            break;
        }
        temp2 = epw;
    }

    /* Short period preliminary quantities */
    let ecose = temp5 + temp6;
    let esine = temp3 - temp4;
    let elsq = axn * axn + ayn * ayn;
    let temp = 1.0 - elsq;
    let pl = a * temp;
    let r = a * (1.0 - ecose);
    let temp1 = 1.0 / r;
    let rdot = XKE * a.sqrt() * esine * temp1;
    let rfdot = XKE * pl.sqrt() * temp1;
    let temp2 = a * temp1;
    let betal = temp.sqrt();
    temp3 = 1.0 / (1.0 + betal);
    let cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
    let sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
    let u = ac_tan(sinu, cosu);
    let sin2u = 2.0 * sinu * cosu;
    let cos2u = 2.0 * cosu * cosu - 1.0;
    let temp = 1.0 / pl;
    let temp1 = CK2 * temp;
    let temp2 = temp1 * temp;

    /* Update for short periodics */
    let rk = r * (1.0 - 1.5 * temp2 * betal * sdp.x3thm1) + 0.5 * temp1 * sdp.x1mth2 * cos2u;
    let uk = u - 0.25 * temp2 * sdp.x7thm1 * sin2u;
    let xnodek = sdp.xnode + 1.5 * temp2 * sdp.cosio * sin2u;
    let xinck = sdp.xinc + 1.5 * temp2 * sdp.cosio * sdp.sinio * cos2u;
    let rdotk = rdot - sdp.xn * temp1 * sdp.x1mth2 * sin2u;
    let rfdotk = rfdot + sdp.xn * temp1 * (sdp.x1mth2 * cos2u + 1.5 * sdp.x3thm1);

    /* Orientation vectors */
    let (sinuk, cosuk) = uk.sin_cos();
    let (sinik, cosik) = xinck.sin_cos();
    let (sinnok, cosnok) = xnodek.sin_cos();
    let xmx = -sinnok * cosik;
    let xmy = cosnok * cosik;
    let ux = xmx * sinuk + cosnok * cosuk;
    let uy = xmy * sinuk + sinnok * cosuk;
    let uz = sinik * sinuk;
    let vx = xmx * cosuk - cosnok * sinuk;
    let vy = xmy * cosuk - sinnok * sinuk;
    let vz = sinik * cosuk;

    /* Position and velocity */
    sdp.pos.0 = rk * ux;
    sdp.pos.1 = rk * uy;
    sdp.pos.2 = rk * uz;
    sdp.vel.0 = rdotk * ux + rfdotk * vx;
    sdp.vel.1 = rdotk * uy + rfdotk * vy;
    sdp.vel.2 = rdotk * uz + rfdotk * vz;

    /* Phase in rads */
    sdp.phase = xlt - sdp.xnode - sdp.omgadf + TWOPI;
    if sdp.phase < 0.0 {
        sdp.phase += TWOPI;
    }
    sdp.phase = fmod2p(sdp.phase);

    tle.omegao1 = sdp.omgadf;
    tle.xincl1 = sdp.xinc;
    tle.xnodeo1 = sdp.xnode;
}

#[derive(Default)]
struct DeepData {
    lunar_terms_done: bool,
    resonance: bool,
    synchronus: bool,
    do_loop: bool,
    epoch_restart: bool,
    thgr: f64,
    xnq: f64,
    xqncl: f64,
    omegaq: f64,
    preep: f64,
    zcosil: f64,
    zsinil: f64,
    zsinhl: f64,
    zcoshl: f64,
    zmol: f64,
    zcosgl: f64,
    zsingl: f64,
    zmos: f64,
    savtsn: f64,
    ee2: f64,
    e3: f64,
    xi2: f64,
    xi3: f64,
    xl2: f64,
    xl3: f64,
    xl4: f64,
    xgh2: f64,
    xgh3: f64,
    xgh4: f64,
    xh2: f64,
    xh3: f64,
    sse: f64,
    ssi: f64,
    ssl: f64,
    ssh: f64,
    ssg: f64,
    se2: f64,
    si2: f64,
    sl2: f64,
    sgh2: f64,
    sh2: f64,
    se3: f64,
    si3: f64,
    sl3: f64,
    sgh3: f64,
    sh3: f64,
    sl4: f64,
    sgh4: f64,
    d2201: f64,
    d2211: f64,
    d3210: f64,
    d3222: f64,
    d4410: f64,
    d4422: f64,
    d5220: f64,
    d5232: f64,
    d5421: f64,
    d5433: f64,
    xlamo: f64,
    del1: f64,
    del2: f64,
    del3: f64,
    fasx2: f64,
    fasx4: f64,
    fasx6: f64,
    xfact: f64,
    xli: f64,
    xni: f64,
    atime: f64,
    stepp: f64,
    stepn: f64,
    step2: f64,
    sghs: f64,
    shs: f64,
    sghl: f64,
    sh1: f64,
    pe: f64,
    pinc: f64,
    pl: f64,
}

/* DEEP */
/* This function is used by SDP4 to add lunar and solar */
/* perturbation effects to deep-space orbit objects.    */
fn deep(ientry: DPEntry, sdp: &mut SDPData, dps: &mut DeepData, tle: &crate::tle::TLE) {
    let [mut a1,mut a2,mut a3,mut a4,mut a5,mut a6,mut a7,mut a8,mut a9,mut a10,ainv2,mut alfdp,aqnv,
         mut sgh,sini2,sinis,sinok,mut sh,mut si,sil,day,mut betdp,dalf,
         mut bfact,c,mut cc,cosis,cosok,cosq,ctem,f322,mut zx,zy,
         dbet,dls,eoc,eq,mut f2,f220,f221,mut f3,f311,f321,xnoh,
         mut f330,f441,f442,f522,f523,f542,f543,g200,g201,
         g211,mut pgh,mut ph,mut s1,mut s2,mut s3,mut s4,mut s5,mut s6,mut s7,mut se,sel,ses,mut xls,
         g300,g310,g322,g410,g422,g520,g521,g532,g533,gam,
         sinq,mut sinzf,sis,mut sl,sll,sls,stem,mut temp,mut temp1,mut x1,mut x2,
         mut x2li,mut x2omi,mut x3,mut x4,mut x5,mut x6,mut x7,mut x8,xl,mut xldot,xmao,mut xnddt,
         mut xndot,xno2,xnodce,xnoi,mut xomi,xpidot,mut z1,mut z11,mut z12,mut z13,
         mut z2,mut z21,mut z22,mut z23,mut z3,mut z31,mut z32,mut z33,mut ze,mut zf,mut zm,mut zn,
         mut zsing,mut zsinh,mut zsini,mut zcosg,mut zcosh,mut zcosi,mut delt,mut ft]: [f64; 130];
    delt = 0.0;
    ft = 0.0;

    match ientry {
        DPEntry::Init => {
            /* Entrance for deep space initialization */
            dps.thgr = theta_g(tle.epoch, &mut sdp.ds50);
            eq = tle.eo;
            dps.xnq = sdp.xnodp;
            aqnv = 1.0 / sdp.aodp;
            dps.xqncl = tle.xincl;
            xmao = tle.xmo;
            xpidot = sdp.omgdot + sdp.xnodot;
            let (sinq_t, cosq_t) = tle.xnodeo.sin_cos();
            sinq = sinq_t;
            cosq = cosq_t;
            dps.omegaq = tle.omegao;
            dps.preep = 0.0;

            /* Initialize lunar solar terms */
            day = sdp.ds50 + 18261.5; /*Days since 1900 Jan 0.5*/
            if day != dps.preep {
                dps.preep = day;
                xnodce = 4.5236020 - 9.2422029E-4 * day;
                let (stem_t, ctem_t) = xnodce.sin_cos();
                stem = stem_t;
                ctem = ctem_t;
                dps.zcosil = 0.91375164 - 0.03568096 * ctem;
                dps.zsinil = (1.0 - dps.zcosil * dps.zcosil).sqrt();
                dps.zsinhl = 0.089683511 * stem / dps.zsinil;
                dps.zcoshl = (1.0 - dps.zsinhl * dps.zsinhl).sqrt();
                c = 4.7199672 + 0.22997150 * day;
                gam = 5.8351514 + 0.0019443680 * day;
                dps.zmol = fmod2p(c - gam);
                zx = 0.39785416 * stem / dps.zsinil;
                zy = dps.zcoshl * ctem + 0.91744867 * dps.zsinhl * stem;
                zx = ac_tan(zx, zy);
                zx = gam + zx - xnodce;
                let (zsingl_t, zcosgl_t) = zx.sin_cos();
                dps.zcosgl = zcosgl_t;
                dps.zsingl = zsingl_t;
                dps.zmos = 6.2565837 + 0.017201977 * day;
                dps.zmos = fmod2p(dps.zmos);
            } /* End if(day != preep) */

            /* Do solar terms */
            dps.savtsn = 1E20;
            zcosg = ZCOSGS;
            zsing = ZSINGS;
            zcosi = ZCOSIS;
            zsini = ZSINIS;
            zcosh = cosq;
            zsinh = sinq;
            cc = C1SS;
            zn = ZNS;
            ze = ZES;
            xnoi = 1.0 / dps.xnq;

            /* Loop breaks when Solar terms are done a second */
            /* time, after Lunar terms are initialized        */
            loop {
                /* Solar terms done again after Lunar terms are done */
                a1 = zcosg * zcosh + zsing * zcosi * zsinh;
                a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
                a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
                a8 = zsing * zsini;
                a9 = zsing * zsinh + zcosg * zcosi * zcosh;
                a10 = zcosg * zsini;
                a2 = sdp.cosio * a7 + sdp.sinio * a8;
                a4 = sdp.cosio * a9 + sdp.sinio * a10;
                a5 = -sdp.sinio * a7 + sdp.cosio * a8;
                a6 = -sdp.sinio * a9 + sdp.cosio * a10;
                x1 = a1 * sdp.cosg + a2 * sdp.sing;
                x2 = a3 * sdp.cosg + a4 * sdp.sing;
                x3 = -a1 * sdp.sing + a2 * sdp.cosg;
                x4 = -a3 * sdp.sing + a4 * sdp.cosg;
                x5 = a5 * sdp.sing;
                x6 = a6 * sdp.sing;
                x7 = a5 * sdp.cosg;
                x8 = a6 * sdp.cosg;
                z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
                z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
                z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
                z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * sdp.eosq;
                z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * sdp.eosq;
                z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * sdp.eosq;
                z11 = -6.0 * a1 * a5 + sdp.eosq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
                z12 = -6.0 * (a1 * a6 + a3 * a5)
                    + sdp.eosq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
                z13 = -6.0 * a3 * a6 + sdp.eosq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
                z21 = 6.0 * a2 * a5 + sdp.eosq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
                z22 = 6.0 * (a4 * a5 + a2 * a6)
                    + sdp.eosq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
                z23 = 6.0 * a4 * a6 + sdp.eosq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
                z1 = z1 + z1 + sdp.betao2 * z31;
                z2 = z2 + z2 + sdp.betao2 * z32;
                z3 = z3 + z3 + sdp.betao2 * z33;
                s3 = cc * xnoi;
                s2 = -0.5 * s3 / sdp.betao;
                s4 = s3 * sdp.betao;
                s1 = -15.0 * eq * s4;
                s5 = x1 * x3 + x2 * x4;
                s6 = x2 * x3 + x1 * x4;
                s7 = x2 * x4 - x1 * x3;
                se = s1 * zn * s5;
                si = s2 * zn * (z11 + z13);
                sl = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * sdp.eosq);
                sgh = s4 * zn * (z31 + z33 - 6.0);
                sh = -zn * s2 * (z21 + z23);
                if dps.xqncl < 5.2359877E-2 {
                    sh = 0.0;
                }
                dps.ee2 = 2.0 * s1 * s6;
                dps.e3 = 2.0 * s1 * s7;
                dps.xi2 = 2.0 * s2 * z12;
                dps.xi3 = 2.0 * s2 * (z13 - z11);
                dps.xl2 = -2.0 * s3 * z2;
                dps.xl3 = -2.0 * s3 * (z3 - z1);
                dps.xl4 = -2.0 * s3 * (-21.0 - 9.0 * sdp.eosq) * ze;
                dps.xgh2 = 2.0 * s4 * z32;
                dps.xgh3 = 2.0 * s4 * (z33 - z31);
                dps.xgh4 = -18.0 * s4 * ze;
                dps.xh2 = -2.0 * s2 * z22;
                dps.xh3 = -2.0 * s2 * (z23 - z21);

                if dps.lunar_terms_done {
                    break;
                }

                /* Do lunar terms */
                dps.sse = se;
                dps.ssi = si;
                dps.ssl = sl;
                dps.ssh = sh / sdp.sinio;
                dps.ssg = sgh - sdp.cosio * dps.ssh;
                dps.se2 = dps.ee2;
                dps.si2 = dps.xi2;
                dps.sl2 = dps.xl2;
                dps.sgh2 = dps.xgh2;
                dps.sh2 = dps.xh2;
                dps.se3 = dps.e3;
                dps.si3 = dps.xi3;
                dps.sl3 = dps.xl3;
                dps.sgh3 = dps.xgh3;
                dps.sh3 = dps.xh3;
                dps.sl4 = dps.xl4;
                dps.sgh4 = dps.xgh4;
                zcosg = dps.zcosgl;
                zsing = dps.zsingl;
                zcosi = dps.zcosil;
                zsini = dps.zsinil;
                zcosh = dps.zcoshl * cosq + dps.zsinhl * sinq;
                zsinh = sinq * dps.zcoshl - cosq * dps.zsinhl;
                zn = ZNL;
                cc = C1L;
                ze = ZEL;
                dps.lunar_terms_done = true;
            } /* End of for(;;) */

            dps.sse = dps.sse + se;
            dps.ssi = dps.ssi + si;
            dps.ssl = dps.ssl + sl;
            dps.ssg = dps.ssg + sgh - sdp.cosio / sdp.sinio * sh;
            dps.ssh = dps.ssh + sh / sdp.sinio;

            /* Geopotential resonance initialization for 12 hour orbits */
            dps.resonance = false;
            dps.synchronus = false;

            if !((dps.xnq < 0.0052359877) && (dps.xnq > 0.0034906585)) {
                if (dps.xnq < 0.00826) || (dps.xnq > 0.00924) {
                    return;
                }
                if eq < 0.5 {
                    return;
                }
                dps.resonance = true;
                eoc = eq * sdp.eosq;
                g201 = -0.306 - (eq - 0.64) * 0.440;
                if eq <= 0.65 {
                    g211 = 3.616 - 13.247 * eq + 16.290 * sdp.eosq;
                    g310 = -19.302 + 117.390 * eq - 228.419 * sdp.eosq + 156.591 * eoc;
                    g322 = -18.9068 + 109.7927 * eq - 214.6334 * sdp.eosq + 146.5816 * eoc;
                    g410 = -41.122 + 242.694 * eq - 471.094 * sdp.eosq + 313.953 * eoc;
                    g422 = -146.407 + 841.880 * eq - 1629.014 * sdp.eosq + 1083.435 * eoc;
                    g520 = -532.114 + 3017.977 * eq - 5740.0 * sdp.eosq + 3708.276 * eoc;
                } else {
                    g211 = -72.099 + 331.819 * eq - 508.738 * sdp.eosq + 266.724 * eoc;
                    g310 = -346.844 + 1582.851 * eq - 2415.925 * sdp.eosq + 1246.113 * eoc;
                    g322 = -342.585 + 1554.908 * eq - 2366.899 * sdp.eosq + 1215.972 * eoc;
                    g410 = -1052.797 + 4758.686 * eq - 7193.992 * sdp.eosq + 3651.957 * eoc;
                    g422 = -3581.69 + 16178.11 * eq - 24462.77 * sdp.eosq + 12422.52 * eoc;
                    if eq <= 0.715 {
                        g520 = 1464.74 - 4664.75 * eq + 3763.64 * sdp.eosq;
                    } else {
                        g520 = -5149.66 + 29936.92 * eq - 54087.36 * sdp.eosq + 31324.56 * eoc;
                    }
                } /* End if (eq <= 0.65) */

                if eq < 0.7 {
                    g533 = -919.2277 + 4988.61 * eq - 9064.77 * sdp.eosq + 5542.21 * eoc;
                    g521 = -822.71072 + 4568.6173 * eq - 8491.4146 * sdp.eosq + 5337.524 * eoc;
                    g532 = -853.666 + 4690.25 * eq - 8624.77 * sdp.eosq + 5341.4 * eoc;
                } else {
                    g533 = -37995.78 + 161616.52 * eq - 229838.2 * sdp.eosq + 109377.94 * eoc;
                    g521 = -51752.104 + 218913.95 * eq - 309468.16 * sdp.eosq + 146349.42 * eoc;
                    g532 = -40023.88 + 170470.89 * eq - 242699.48 * sdp.eosq + 115605.82 * eoc;
                } /* End if (eq <= 0.7) */

                sini2 = sdp.sinio * sdp.sinio;
                f220 = 0.75 * (1.0 + 2.0 * sdp.cosio + sdp.theta2);
                f221 = 1.5 * sini2;
                f321 = 1.875 * sdp.sinio * (1.0 - 2.0 * sdp.cosio - 3.0 * sdp.theta2);
                f322 = -1.875 * sdp.sinio * (1.0 + 2.0 * sdp.cosio - 3.0 * sdp.theta2);
                f441 = 35.0 * sini2 * f220;
                f442 = 39.3750 * sini2 * sini2;
                f522 = 9.84375
                    * sdp.sinio
                    * (sini2 * (1.0 - 2.0 * sdp.cosio - 5.0 * sdp.theta2)
                        + 0.33333333 * (-2.0 + 4.0 * sdp.cosio + 6.0 * sdp.theta2));
                f523 = sdp.sinio
                    * (4.92187512 * sini2 * (-2.0 - 4.0 * sdp.cosio + 10.0 * sdp.theta2)
                        + 6.56250012 * (1.0 + 2.0 * sdp.cosio - 3.0 * sdp.theta2));
                f542 = 29.53125
                    * sdp.sinio
                    * (2.0 - 8.0 * sdp.cosio
                        + sdp.theta2 * (-12.0 + 8.0 * sdp.cosio + 10.0 * sdp.theta2));
                f543 = 29.53125
                    * sdp.sinio
                    * (-2.0 - 8.0 * sdp.cosio
                        + sdp.theta2 * (12.0 + 8.0 * sdp.cosio - 10.0 * sdp.theta2));
                xno2 = dps.xnq * dps.xnq;
                ainv2 = aqnv * aqnv;
                temp1 = 3.0 * xno2 * ainv2;
                temp = temp1 * ROOT22;
                dps.d2201 = temp * f220 * g201;
                dps.d2211 = temp * f221 * g211;
                temp1 = temp1 * aqnv;
                temp = temp1 * ROOT32;
                dps.d3210 = temp * f321 * g310;
                dps.d3222 = temp * f322 * g322;
                temp1 = temp1 * aqnv;
                temp = 2.0 * temp1 * ROOT44;
                dps.d4410 = temp * f441 * g410;
                dps.d4422 = temp * f442 * g422;
                temp1 = temp1 * aqnv;
                temp = temp1 * ROOT52;
                dps.d5220 = temp * f522 * g520;
                dps.d5232 = temp * f523 * g532;
                temp = 2.0 * temp1 * ROOT54;
                dps.d5421 = temp * f542 * g521;
                dps.d5433 = temp * f543 * g533;
                dps.xlamo = xmao + tle.xnodeo + tle.xnodeo - dps.thgr - dps.thgr;
                bfact = sdp.xmdot + sdp.xnodot + sdp.xnodot - THDT - THDT;
                bfact = bfact + dps.ssl + dps.ssh + dps.ssh;
            } else {
                dps.resonance = true;
                dps.synchronus = true;
                /* Synchronous resonance terms initialization */
                g200 = 1.0 + sdp.eosq * (-2.5 + 0.8125 * sdp.eosq);
                g310 = 1.0 + 2.0 * sdp.eosq;
                g300 = 1.0 + sdp.eosq * (-6.0 + 6.60937 * sdp.eosq);
                f220 = 0.75 * (1.0 + sdp.cosio) * (1.0 + sdp.cosio);
                f311 = 0.9375 * sdp.sinio * sdp.sinio * (1.0 + 3.0 * sdp.cosio)
                    - 0.75 * (1.0 + sdp.cosio);
                f330 = 1.0 + sdp.cosio;
                f330 = 1.875 * f330 * f330 * f330;
                dps.del1 = 3.0 * dps.xnq * dps.xnq * aqnv * aqnv;
                dps.del2 = 2.0 * dps.del1 * f220 * g200 * Q22;
                dps.del3 = 3.0 * dps.del1 * f330 * g300 * Q33 * aqnv;
                dps.del1 = dps.del1 * f311 * g310 * Q31 * aqnv;
                dps.fasx2 = 0.13130908;
                dps.fasx4 = 2.8843198;
                dps.fasx6 = 0.37448087;
                dps.xlamo = xmao + tle.xnodeo + tle.omegao - dps.thgr;
                bfact = sdp.xmdot + xpidot - THDT;
                bfact = bfact + dps.ssl + dps.ssg + dps.ssh;
            }

            dps.xfact = bfact - dps.xnq;

            /* Initialize integrator */
            dps.xli = dps.xlamo;
            dps.xni = dps.xnq;
            dps.atime = 0.0;
            dps.stepp = 720.0;
            dps.stepn = -720.0;
            dps.step2 = 259200.0;
            /* End case dpinit: */
            return;
        }

        DPEntry::Sec => {
            /* Entrance for deep space secular effects */
            sdp.xll = sdp.xll + dps.ssl * sdp.t;
            sdp.omgadf = sdp.omgadf + dps.ssg * sdp.t;
            sdp.xnode = sdp.xnode + dps.ssh * sdp.t;
            sdp.em = tle.eo + dps.sse * sdp.t;
            sdp.xinc = tle.xincl + dps.ssi * sdp.t;
            if sdp.xinc < 0.0 {
                sdp.xinc = -sdp.xinc;
                sdp.xnode = sdp.xnode + PI;
                sdp.omgadf = sdp.omgadf - PI;
            }
            if !dps.resonance {
                return;
            }

            loop {
                if (dps.atime == 0.0)
                    || ((sdp.t >= 0.0) && (dps.atime < 0.0))
                    || ((sdp.t < 0.0) && (dps.atime >= 0.0))
                {
                    /* Epoch restart */
                    if sdp.t >= 0.0 {
                        delt = dps.stepp;
                    } else {
                        delt = dps.stepn;
                    }

                    dps.atime = 0.0;
                    dps.xni = dps.xnq;
                    dps.xli = dps.xlamo;
                } else {
                    if sdp.t.abs() >= dps.atime.abs() {
                        if sdp.t > 0.0 {
                            delt = dps.stepp;
                        } else {
                            delt = dps.stepn;
                        }
                    }
                }

                loop {
                    if (sdp.t - dps.atime).abs() >= dps.stepp {
                        dps.do_loop = true;
                        dps.epoch_restart = false;
                    } else {
                        ft = sdp.t - dps.atime;
                        dps.do_loop = false;
                    }

                    if sdp.t.abs() < dps.atime.abs() {
                        if sdp.t >= 0.0 {
                            delt = dps.stepn;
                        } else {
                            delt = dps.stepp;
                        }
                        dps.do_loop = true;
                        dps.epoch_restart = true;
                    }

                    /* Dot terms calculated */
                    if dps.synchronus {
                        xndot = dps.del1 * (dps.xli - dps.fasx2).sin()
                            + dps.del2 * (2.0 * (dps.xli - dps.fasx4)).sin()
                            + dps.del3 * (3.0 * (dps.xli - dps.fasx6)).sin();
                        xnddt = dps.del1 * (dps.xli - dps.fasx2).cos()
                            + 2.0 * dps.del2 * (2.0 * (dps.xli - dps.fasx4)).cos()
                            + 3.0 * dps.del3 * (3.0 * (dps.xli - dps.fasx6)).cos();
                    } else {
                        xomi = dps.omegaq + sdp.omgdot * dps.atime;
                        x2omi = xomi + xomi;
                        x2li = dps.xli + dps.xli;
                        xndot = dps.d2201 * (x2omi + dps.xli - G22).sin()
                            + dps.d2211 * (dps.xli - G22).sin()
                            + dps.d3210 * (xomi + dps.xli - G32).sin()
                            + dps.d3222 * (-xomi + dps.xli - G32).sin()
                            + dps.d4410 * (x2omi + x2li - G44).sin()
                            + dps.d4422 * (x2li - G44).sin()
                            + dps.d5220 * (xomi + dps.xli - G52).sin()
                            + dps.d5232 * (-xomi + dps.xli - G52).sin()
                            + dps.d5421 * (xomi + x2li - G54).sin()
                            + dps.d5433 * (-xomi + x2li - G54).sin();
                        xnddt = dps.d2201 * (x2omi + dps.xli - G22).cos()
                            + dps.d2211 * (dps.xli - G22).cos()
                            + dps.d3210 * (xomi + dps.xli - G32).cos()
                            + dps.d3222 * (-xomi + dps.xli - G32).cos()
                            + dps.d5220 * (xomi + dps.xli - G52).cos()
                            + dps.d5232 * (-xomi + dps.xli - G52).cos()
                            + 2.0
                                * (dps.d4410 * (x2omi + x2li - G44).cos()
                                    + dps.d4422 * (x2li - G44).cos()
                                    + dps.d5421 * (xomi + x2li - G54).cos()
                                    + dps.d5433 * (-xomi + x2li - G54))
                                    .cos();
                    } /* End of if (isFlagSet(SYNCHRONOUS_FLAG)) */

                    xldot = dps.xni + dps.xfact;
                    xnddt = xnddt * xldot;

                    if dps.do_loop {
                        dps.xli = dps.xli + xldot * delt + xndot * dps.step2;
                        dps.xni = dps.xni + xndot * delt + xnddt * dps.step2;
                        dps.atime = dps.atime + delt;
                    }

                    if !((dps.do_loop) && (!dps.epoch_restart)) {
                        break;
                    }
                }

                if !((dps.do_loop) && (dps.epoch_restart)) {
                    break;
                }
            }

            sdp.xn = dps.xni + xndot * ft + xnddt * ft * ft * 0.5;
            xl = dps.xli + xldot * ft + xndot * ft * ft * 0.5;
            temp = -sdp.xnode + dps.thgr + sdp.t * THDT;

            if dps.synchronus {
                sdp.xll = xl + temp + temp;
            } else {
                sdp.xll = xl - sdp.omgadf + temp;
            }

            return;
            /*End case dpsec: */
        }

        DPEntry::Per => {
            /* Entrance for lunar-solar periodics */
            let (sinis_t, cosis_t) = sdp.xinc.sin_cos();
            sinis = sinis_t;
            cosis = cosis_t;
            if (dps.savtsn - sdp.t).abs() >= 30.0 {
                dps.savtsn = sdp.t;
                zm = dps.zmos + ZNS * sdp.t;
                zf = zm + 2.0 * ZES * zm.sin();
                sinzf = zf.sin();
                f2 = 0.5 * sinzf * sinzf - 0.25;
                f3 = -0.5 * sinzf * zf.cos();
                ses = dps.se2 * f2 + dps.se3 * f3;
                sis = dps.si2 * f2 + dps.si3 * f3;
                sls = dps.sl2 * f2 + dps.sl3 * f3 + dps.sl4 * sinzf;
                dps.sghs = dps.sgh2 * f2 + dps.sgh3 * f3 + dps.sgh4 * sinzf;
                dps.shs = dps.sh2 * f2 + dps.sh3 * f3;
                zm = dps.zmol + ZNL * sdp.t;
                zf = zm + 2.0 * ZEL * zm.sin();
                sinzf = zf.sin();
                f2 = 0.5 * sinzf * sinzf - 0.25;
                f3 = -0.5 * sinzf * zf.cos();
                sel = dps.ee2 * f2 + dps.e3 * f3;
                sil = dps.xi2 * f2 + dps.xi3 * f3;
                sll = dps.xl2 * f2 + dps.xl3 * f3 + dps.xl4 * sinzf;
                dps.sghl = dps.xgh2 * f2 + dps.xgh3 * f3 + dps.xgh4 * sinzf;
                dps.sh1 = dps.xh2 * f2 + dps.xh3 * f3;
                dps.pe = ses + sel;
                dps.pinc = sis + sil;
                dps.pl = sls + sll;
            }

            pgh = dps.sghs + dps.sghl;
            ph = dps.shs + dps.sh1;
            sdp.xinc = sdp.xinc + dps.pinc;
            sdp.em = sdp.em + dps.pe;

            if dps.xqncl >= 0.2 {
                /* Apply periodics directly */
                ph = ph / sdp.sinio;
                pgh = pgh - sdp.cosio * ph;
                sdp.omgadf = sdp.omgadf + pgh;
                sdp.xnode = sdp.xnode + ph;
                sdp.xll = sdp.xll + dps.pl;
            } else {
                /* Apply periodics with Lyddane modification */
                let (sinok_t, cosok_t) = sdp.xnode.sin_cos();
                sinok = sinok_t;
                cosok = cosok_t;
                alfdp = sinis * sinok;
                betdp = sinis * cosok;
                dalf = ph * cosok + dps.pinc * cosis * sinok;
                dbet = -ph * sinok + dps.pinc * cosis * cosok;
                alfdp = alfdp + dalf;
                betdp = betdp + dbet;
                sdp.xnode = fmod2p(sdp.xnode);
                xls = sdp.xll + sdp.omgadf + cosis * sdp.xnode;
                dls = dps.pl + pgh - dps.pinc * sdp.xnode * sinis;
                xls = xls + dls;
                xnoh = sdp.xnode;
                sdp.xnode = ac_tan(alfdp, betdp);

                /* This is a patch to Lyddane modification */
                /* suggested by Rob Matson. */
                if (xnoh - sdp.xnode).abs() > PI {
                    if sdp.xnode < xnoh {
                        sdp.xnode += TWOPI;
                    } else {
                        sdp.xnode -= TWOPI;
                    }
                }

                sdp.xll = sdp.xll + dps.pl;
                sdp.omgadf = xls - sdp.xll - sdp.xinc.cos() * sdp.xnode;
            }
            return;
        }
    }
}
