use super::*;

const TLE1: &str = "TEST SAT SGP 001
1 88888U 88888A   80275.98708465  .00073094  13844-3  66816-4 0  5559
2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   103";

const TLE2: &str = "TEST SAT SDP 001
1 11801U 88888A   80230.29629788  .01431103  00000-0  14311-1 0  5552
2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848   102";

struct DataSet {
    t: f64,
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

const DATA1: [DataSet; 5] = [
    DataSet {
        t: 0.0,
        x: 2328.97048951,
        y: -5995.22076416,
        z: 1719.97067261,
        vx: 2.91207230,
        vy: -0.98341546,
        vz: -7.09081703,
    },
    DataSet {
        t: 360.0,
        x: 2456.10705566,
        y: -6071.93853760,
        z: 1222.89727783,
        vx: 2.67938992,
        vy: -0.44829041,
        vz: -7.22879231,
    },
    DataSet {
        t: 720.0,
        x: 2567.56195068,
        y: -6112.50384522,
        z: 713.96397400,
        vx: 2.44024599,
        vy: 0.09810869,
        vz: -7.31995916,
    },
    DataSet {
        t: 1080.0,
        x: 2663.09078980,
        y: -6115.48229980,
        z: 196.39640427,
        vx: 2.19611958,
        vy: 0.65241995,
        vz: -7.36282432,
    },
    DataSet {
        t: 1440.0,
        x: 2742.55133057,
        y: -6079.67144775,
        z: -326.38095856,
        vx: 1.94850229,
        vy: 1.21106251,
        vz: -7.35619372,
    },
];

const DATA2: [DataSet; 5] = [
    DataSet {
        t: 0.0,
        x: 7473.37066650,
        y: 428.95261765,
        z: 5828.74786377,
        vx: 5.1071513,
        vy: 6.44468284,
        vz: -0.18613096,
    },
    DataSet {
        t: 360.0,
        x: -3305.22537232,
        y: 32410.86328125,
        z: -24697.17675781,
        vx: -1.30113538,
        vy: -1.15131518,
        vz: -0.28333528,
    },
    DataSet {
        t: 720.0,
        x: 14271.28759766,
        y: 24110.46411133,
        z: -4725.76837158,
        vx: -0.32050445,
        vy: 2.67984074,
        vz: -2.08405289,
    },
    DataSet {
        t: 1080.0,
        x: -9990.05883789,
        y: 22717.35522461,
        z: -23616.890662501,
        vx: -1.01667246,
        vy: -2.29026759,
        vz: 0.72892364,
    },
    DataSet {
        t: 1440.0,
        x: 9787.86975097,
        y: 33753.34667969,
        z: -15030.81176758,
        vx: -1.09425966,
        vy: 0.92358845,
        vz: -1.52230928,
    },
];
#[test]
fn test_sgp() {
    let mut tle = crate::tle::TLE::parse(TLE1);
    let mut sat = Sat::new();

    for DataSet {
        t,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    } in &DATA1
    {
        sat.update(&mut tle, *t);
        let crate::coords::EarthCenteredInertial {
            x: c_x,
            y: c_y,
            z: c_z,
        } = sat.get_coords().unwrap();
        let crate::coords::eci::ECIVector {
            x: c_vx,
            y: c_vy,
            z: c_vz,
        } = sat.get_vel().unwrap();

        assert_approx_eq!(x, c_x, 0.5);
        assert_approx_eq!(y, c_y, 0.5);
        assert_approx_eq!(z, c_z, 0.5);
        assert_approx_eq!(vx, c_vx, 0.5);
        assert_approx_eq!(vy, c_vy, 0.5);
        assert_approx_eq!(vz, c_vz, 0.5);
    }
}

#[test]
fn test_sdp() {
    let mut tle = crate::tle::TLE::parse(TLE2);
    let mut sat = Sat::new();

    for DataSet {
        t,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    } in &DATA2
    {
        sat.update(&mut tle, *t);
        let crate::coords::EarthCenteredInertial {
            x: c_x,
            y: c_y,
            z: c_z,
        } = sat.get_coords().unwrap();
        let crate::coords::eci::ECIVector {
            x: c_vx,
            y: c_vy,
            z: c_vz,
        } = sat.get_vel().unwrap();

        assert_approx_eq!(x, c_x, 0.5);
        assert_approx_eq!(y, c_y, 0.5);
        assert_approx_eq!(z, c_z, 0.5);
        assert_approx_eq!(vx, c_vx, 0.5);
        assert_approx_eq!(vy, c_vy, 0.5);
        assert_approx_eq!(vz, c_vz, 0.5);
    }
}
