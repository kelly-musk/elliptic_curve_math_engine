//!## Modular Arithmetic
//!### Implement Add, sub, mul and div for secp256k1 U256: ([[u64;4]]) data type

use primitive_types::{U256, U512};
use std::ops::{Add, Div, Mul, Sub};


/// Prime of the secp256k1 curve
///
/// U256: ([[u64;4]])
///
/// y^2 = x^3 + 7
pub const P: U256 = U256([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);

/// FieldElement which would be the basis of our curve points
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct FieldElement {
    pub value: U256,
}

impl FieldElement {
    /// Create new instance Of FieldElement type
    pub(crate) fn new(value: U256) -> Self {
        let mut res = value % P;
        if res < U256::zero() {
            res += P;
        }
        Self { value: res }
    }

    /// Using Extended Euclidean Algorithm: `ax + by = gcd(a,b)` to find inverse
    ///
    /// If gcd(a,m) == 1, then ax + my = 1, so x is the modular inverse of a mod m / Prime field
    pub(crate) fn inverse(&self) -> Self {
        // To ensure the value is not zero 
        if self.value == U256::zero() {
            panic!("Cannot inverse a zero value");
        }
        // it should be from the field i.e (1 ..= P-1)
        if self.value >= P {
            panic!("the value {:?} is not in the field {:?}", self.value, P);
        }

        // Extended Euclidean Algorithm with unsigned arithmetic
        let (mut t, mut new_t) = (U256::zero(), U256::one());
        let (mut r, mut new_r) = (P, self.value);

        while new_r != U256::zero() {
            let quotient = r / new_r;

            // Update t: handle subtraction that might go negative
            // Instead of t - quotient * new_t, we compute it modulo P
            let prod = multiply(quotient, new_t);
            let next_t = if t >= prod {
                t - prod
            } else {
                P - (prod - t)
            };
            (t, new_t) = (new_t, next_t);
            (r, new_r) = (new_r, r - multiply(quotient, new_r));
        }

        if r > U256::one() {
            panic!("the inverse does not exist");
        }

        FieldElement::new(t)
    }
}


/// Helper function to handle multiplication for U256 values and avoid overflows
/// 
/// This converts to a u512 which even maxU256 ^ 2 can never overflow, then performs modulo of P in u512 form, 
/// Then takes the least sig bits i.e. little endian and converts back to a U256([[u64;4]])
pub(crate) fn multiply(a: U256, b: U256) -> U256 {
    let result = a.full_mul(b);
    let reduced = result % U512::from(P);
    let lower_256 = U256([reduced.0[0], reduced.0[1], reduced.0[2], reduced.0[3]]);
    lower_256
}

// Set various arithmetic for the field points 
impl Add for FieldElement {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        FieldElement::new(self.value + other.value)
    }
}

impl Sub for FieldElement {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let res = if self.value >= other.value {
            self.value - other.value
        } else {
            P - (other.value - self.value)
        };
        FieldElement::new(res)
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        // Use full_mul to handle 256-bit * 256-bit = 512-bit multiplication
        let result = self.value.full_mul(other.value);
        // Reduce the 512-bit result modulo P
        let reduced = result % primitive_types::U512::from(P);
        // Extract lower 256 bits from U512
        let lower_256 = U256([reduced.0[0], reduced.0[1], reduced.0[2], reduced.0[3]]);
        FieldElement::new(lower_256)
    }
}

impl Div for FieldElement {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        FieldElement::new(multiply(self.value, other.inverse().value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_element_creation() {
        let a = FieldElement::new(U256::from(5));
        assert_eq!(a.value, U256::from(5));
    }

    #[test]
    fn test_field_element_modular_reduction() {
        // Test that values >= P are reduced modulo P
        let a = FieldElement::new(P + U256::from(10));
        assert_eq!(a.value, U256::from(10));
    }

    #[test]
    fn test_add_simple() {
        let a = FieldElement::new(U256::from(5));
        let b = FieldElement::new(U256::from(7));
        let result = a + b;
        assert_eq!(result.value, U256::from(12));
    }

    #[test]
    fn test_add_with_modular_wrap() {
        // Test addition that wraps around the modulus
        let a = FieldElement::new(P - U256::from(5));
        let b = FieldElement::new(U256::from(10));
        let result = a + b;
        assert_eq!(result.value, U256::from(5));
    }

    #[test]
    fn test_sub_simple() {
        let a = FieldElement::new(U256::from(10));
        let b = FieldElement::new(U256::from(3));
        let result = a - b;
        assert_eq!(result.value, U256::from(7));
    }

    #[test]
    fn test_sub_with_modular_wrap() {
        // Test subtraction that wraps around (negative result)
        let a = FieldElement::new(U256::from(5));
        let b = FieldElement::new(U256::from(10));
        let result = a - b;
        // Result should be P - 5
        assert_eq!(result.value, P - U256::from(5));
    }

    #[test]
    fn test_mul_simple() {
        let a = FieldElement::new(U256::from(6));
        let b = FieldElement::new(U256::from(7));
        let result = a * b;
        assert_eq!(result.value, U256::from(42));
    }

    #[test]
    fn test_mul_with_modular_reduction() {
        // Test multiplication that requires modular reduction
        let a = FieldElement::new(P - U256::from(1));
        let b = FieldElement::new(U256::from(2));
        let result = a * b;
        // (P - 1) * 2 = 2P - 2 ≡ P - 2 (mod P)
        assert_eq!(result.value, P - U256::from(2));
    }

    #[test]
    fn test_inverse_simple() {
        // Test that 2 * inverse(2) ≡ 1 (mod P)
        let a = FieldElement::new(U256::from(2));
        let a_inv = a.inverse();
        let result = FieldElement::new(U256::from(2)) * a_inv;
        assert_eq!(result.value, U256::from(1));
    }

    #[test]
    fn test_inverse_larger_value() {
        // Test inverse of a larger value
        let a = FieldElement::new(U256::from(12345));
        let a_inv = a.inverse();
        let result = FieldElement::new(U256::from(12345)) * a_inv;
        assert_eq!(result.value, U256::from(1));
    }

    #[test]
    #[should_panic(expected = "Cannot inverse a zero value")]
    fn test_inverse_zero_panics() {
        let zero = FieldElement::new(U256::zero());
        zero.inverse();
    }

    #[test]
    fn test_div_simple() {
        // Test that 10 / 2 = 5
        let a = FieldElement::new(U256::from(10));
        let b = FieldElement::new(U256::from(2));
        let result = a / b;
        assert_eq!(result.value, U256::from(5));
    }

    #[test]
    fn test_div_with_inverse() {
        // Test that (a / b) * b = a
        let a = FieldElement::new(U256::from(42));
        let b = FieldElement::new(U256::from(7));
        let result = a / b;
        let check = result * FieldElement::new(U256::from(7));
        assert_eq!(check.value, U256::from(42));
    }

    #[test]
    #[should_panic(expected = "Cannot inverse a zero value")]
    fn test_div_by_zero_panics() {
        let a = FieldElement::new(U256::from(10));
        let zero = FieldElement::new(U256::zero());
        let _result = a / zero;
    }

    #[test]
    fn test_additive_identity() {
        // Test that a + 0 = a
        let a = FieldElement::new(U256::from(123));
        let zero = FieldElement::new(U256::zero());
        let result = a + zero;
        assert_eq!(result.value, U256::from(123));
    }

    #[test]
    fn test_multiplicative_identity() {
        // Test that a * 1 = a
        let a = FieldElement::new(U256::from(123));
        let one = FieldElement::new(U256::one());
        let result = a * one;
        assert_eq!(result.value, U256::from(123));
    }
}
