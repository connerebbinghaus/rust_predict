/** Table of constant values **/
pub const DE2RA: f64 = 1.74532925E-2; /* Degrees to Radians */
pub const PI: f64 = 3.1415926535898; /* Pi */
pub const PIO2: f64 = 1.5707963267949; /* Pi/2 */
pub const X3PIO2: f64 = 4.71238898; /* 3*Pi/2 */
pub const TWOPI: f64 = 6.2831853071796; /* 2*Pi  */
pub const E6A: f64 = 1.0E-6;
pub const TOTHRD: f64 = 6.6666667E-1; /* 2/3 */
// pub const XJ2: f64 =       1.0826158E-3;   /* J2 Harmonic */
pub const XJ3: f64 = -2.53881E-6; /* J3 Harmonic */
// pub const XJ4: f64 =      -1.65597E-6;     /* J4 Harmonic */
pub const XKE: f64 = 7.43669161E-2;
pub const XKMPER: f64 = 6.378135E3; /* Earth radius km */
pub const XMNPDA: f64 = 1.44E3; /* Minutes per day */
pub const AE: f64 = 1.0;
pub const CK2: f64 = 5.413079E-4;
pub const CK4: f64 = 6.209887E-7;
pub const __F: f64 = 3.352779E-3;
// pub const GE: f64 =        3.986008E5;
pub const __S__: f64 = 1.012229;
pub const QOMS2T: f64 = 1.880279E-09;
pub const SECDAY: f64 = 8.6400E4; /* Seconds per day */
// pub const OMEGA_E: f64 =   1.0027379;
// pub const OMEGA_ER: f64 =  6.3003879;
pub const ZNS: f64 = 1.19459E-5;
pub const C1SS: f64 = 2.9864797E-6;
pub const ZES: f64 = 1.675E-2;
pub const ZNL: f64 = 1.5835218E-4;
pub const C1L: f64 = 4.7968065E-7;
pub const ZEL: f64 = 5.490E-2;
pub const ZCOSIS: f64 = 9.1744867E-1;
pub const ZSINIS: f64 = 3.9785416E-1;
pub const ZSINGS: f64 = -9.8088458E-1;
pub const ZCOSGS: f64 = 1.945905E-1;
// pub const ZCOSHS: f64 =    1.0;
// pub const ZSINHS: f64 =    0.0;
pub const Q22: f64 = 1.7891679E-6;
pub const Q31: f64 = 2.1460748E-6;
pub const Q33: f64 = 2.2123015E-7;
pub const G22: f64 = 5.7686396;
pub const G32: f64 = 9.5240898E-1;
pub const G44: f64 = 1.8014998;
pub const G52: f64 = 1.0508330;
pub const G54: f64 = 4.4108898;
pub const ROOT22: f64 = 1.7891679E-6;
pub const ROOT32: f64 = 3.7393792E-7;
pub const ROOT44: f64 = 7.3636953E-9;
pub const ROOT52: f64 = 1.1428639E-7;
pub const ROOT54: f64 = 2.1765803E-9;
pub const THDT: f64 = 4.3752691E-3;
// pub const RHO: f64 =       1.5696615E-1;
// pub const MFACTOR: f64 =   7.292115E-5;
pub const __SR__: f64 = 6.96000E5; /*Solar radius - kilometers (IAU 76) */
// pub const AU: f64 =        1.49597870E8;   /*Astronomical unit - kilometers (IAU 76) */
/* Entry points of Deep() */
pub enum DPEntry {
    Init,
    Sec,
    Per,
}
