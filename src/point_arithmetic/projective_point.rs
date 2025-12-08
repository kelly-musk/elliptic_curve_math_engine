//! Projective Point (3d Ecpoint including flag(z) for infinity)
//! 
//! This would be used for calculation due to the avoidance of the inverse / division cost
//! 
//! EcPoint for logging and display , ProjectivePoint for calculation

use primitive_types::U256;

use crate::point_arithmetic::{
    modular_arithmetic::FieldElement, 
    point::EcPoint
};

/// Projective Point (X, Y, Z)
/// Represents (X/Z, Y/Z) in Affine coordinates
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectivePoint {
    pub x: FieldElement,
    pub y: FieldElement,
    pub z: FieldElement,
}

impl ProjectivePoint{
    /// To ensure that the point is at infinity, z should be zero
    pub fn is_infinity(&self) -> bool {
        if self.z.value == U256::zero() {
            return true;
        }
        false
    }

    /// When setting to infinity, z == 0 , x == 0 , y can be anything
    /// 
    /// Y^2.Z = X^3.Z + 7.Z^3
    pub fn infinity() -> Self {
        Self {
            x: FieldElement::new(U256::zero()),
            y: FieldElement::new(U256::from(1)),
            z: FieldElement::new(U256::zero()),
        }
    }
}

impl From<EcPoint> for ProjectivePoint {
    /// Convert from EcPoint to Projective
    /// Formula: (x, y) -> (x, y, 1)
    /// Infinity -> (0, 1, 0)
    fn from(ep: EcPoint) -> Self {
        match ep {
            EcPoint::Infinity => Self::infinity(),
            EcPoint::Point { x, y } => Self {
                x,
                y,
                z: FieldElement::new(U256::from(1)),
            },
        }
    }
}

impl From<ProjectivePoint> for EcPoint {
    /// Convert from Projective to Affine
    /// Formula: (X, Y, Z) -> (X/Z, Y/Z)
    /// If Z == 0, return Infinity
    fn from(pp: ProjectivePoint) -> Self {
        // if the z is zero , simply return inifinity here
        if pp.is_infinity() {
            return EcPoint::Infinity;
        }
        // find inverse of z
        let z_inv = pp.z.inverse();
        // return the point (x/z, y/z)
        EcPoint::Point {
            x: pp.x * z_inv,
            y: pp.y * z_inv,
        }
    }
}