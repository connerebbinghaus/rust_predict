pub mod azel;
pub mod eci;
pub mod latlong;
pub mod th;

pub use eci::EarthCenteredInertial;

pub fn calculate_look(
    obj: EarthCenteredInertial,
    obs: latlong::LatitudeLongitude,
    date: chrono::DateTime<chrono::Utc>,
) -> azel::AzElVector {
    let obs_eci = obs.to_eci(date);
    let vec = obj - obs_eci;
    let th = vec.to_topo_horiz(&obs, date);
    th.to_azel()
}
