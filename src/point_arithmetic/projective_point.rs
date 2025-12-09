//! Projective Point (3d Ecpoint including flag(z) for infinity)
//! 
//! This would be used for calculation due to the avoidance of the inverse / division cost
//! 
//! EcPoint for logging and display , JacobianPoint for calculation
//! 
//! Using Jacobian co-ordinates (X, Y, Z) to represent (X/Z^2, Y/Z^3) in EcPoint(x,y) coordinates

use primitive_types::U256;

use crate::point_arithmetic::{
    A, B, modular_arithmetic::FieldElement, multiply, point::EcPoint
};

/// Projective Point (X, Y, Z)
/// Represents (X/Z^2, Y/Z^3) in Affine coordinates
/// 
/// Y^2 = X^3 + aXZ^4 + bZ^6
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JacobianPoint {
    pub x: FieldElement,
    pub y: FieldElement,
    pub z: FieldElement,
}

impl JacobianPoint{
    /// To ensure that the point is at infinity, z should be zero
    pub(crate) fn is_infinity(&self) -> bool {
        if self.z.value == U256::zero() {
            return true;
        }
        false
    }

    /// When setting to infinity, z == 0 , x == 0 , y can be anything
    /// 
    /// Y^2 = X^3 + aXZ^4 + bZ^6
    /// 
    /// Y^2.Z = X^3.Z + 7.Z^3
    pub(crate) fn infinity() -> Self {
        Self {
            x: FieldElement::new(U256::zero()),
            y: FieldElement::new(U256::from(1)),
            z: FieldElement::new(U256::zero()),
        }
    }

    /// Addition for jacobianPoint
    pub(crate) fn add(&self, other: &Self) -> Self{
        if self.is_infinity() {
            return *other;
        }
        if other.is_infinity() {
            return *self;
        }
        if self == other {
            return self.double();
        }
        // sub x = X/Z^2 && y = Y/Z^3
        let z1_square = self.z * self.z;
        let z2_square = other.z * other.z;
        // u1 = X1.Z2^2 
        let u1 = self.x * z2_square;
        // u2 = X2.Z1^2
        let u2 = other.x * z1_square;
        // s1 = Y1.Z2^3
        let s1 = self.y * z2_square * other.z;
        // s2 = Y2.Z1^3
        let s2 = other.y * z1_square * self.z;
        // h = u2 - u1 (change in x)
        let h = u2 - u1;
        // r = s2 - s1 (change in y)
        let r = s2 - s1;
        // x3 = r^2 - h^3 - 2.u1.h^2
        let x3 = (r * r) - (h * h * h) - (FieldElement::new(U256::from(2)) * u1 * h * h);
        // y3 = r.(u1.h^2 - x3) - s1.h^3
        let y3 = r * (u1 * (h * h) - x3) - s1 * (h * h * h);
        // z3 = h.z1.z2
        let z3 = h * self.z * other.z;
        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub(crate) fn double(&self) -> Self {
        // calculate intermediate values / squares
        let a = self.x * self.x;
        let b = self.y * self.y;
        let c = b * self.y;
        // calculate the slope
        let s = FieldElement::new(U256::from(2)) * (((self.x + b) * (self.x + b)) - a - c);
        // calculate the slope numerator
        let m = FieldElement::new(U256::from(3)) * a ;
        let x3 = (m * m) - (FieldElement::new(U256::from(2)) * s);
        let y3 = m * (s - self.x) - (FieldElement::new(U256::from(8)) * c);
        let z3 = FieldElement::new(U256::from(2)) * self.y * self.z;
        // return the new point
        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub(crate) fn scalar_mul(&self, scalar: U256) -> Self {
        todo!()
    }

    pub(crate) fn scalar_div(&self, scalar: U256) -> Self {
        todo!()
    }

    pub(crate) fn sub(&self, other: &Self) -> Self {
        todo!()
    }
}

impl From<EcPoint> for JacobianPoint {
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

impl From<JacobianPoint> for EcPoint {
    /// Convert from Projective to Affine
    /// Formula: (X, Y, Z) -> (X/Z^2, Y/Z^3)
    /// If Z == 0, return Infinity
    fn from(jp: JacobianPoint) -> Self {
        // if the z is zero , simply return inifinity here
        if jp.is_infinity() {
            return EcPoint::Infinity;
        }
        // Find z^2 and z^3
        let z_squared = multiply(jp.z.value, jp.z.value);
        let z_cubed = multiply(z_squared, jp.z.value);
        // find the inverse of both
        let z_squared_inv = FieldElement::new(U256::from(z_squared));
        let z_cubed_inv = FieldElement::new(U256::from(z_cubed));
        // return the point (x/z^2, y/z^3)
        EcPoint::Point {
            x: jp.x * z_squared_inv,
            y: jp.y * z_cubed_inv,
        }
    }
}