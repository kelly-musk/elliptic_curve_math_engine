//!## Point Arithmetic
//!### Implement Point Arithmetic for secp256k1

use primitive_types::U256;
use hex_literal;

use super::{FieldElement, multiply};

/// The weierstrass formula used here is `y^2 = x^3 + 7`
///
/// Original weierstrass formula is `y^2 = x^3 + ax + b`
///
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Point {
    x: FieldElement,
    y: FieldElement,
}
pub const A: U256 = U256::zero();
pub const B: U256 = U256([7, 0, 0, 0]);
/// Gx coordiante for Generator point
pub const G_X_BYTES: [u8;32] = hex_literal::hex!("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798");
/// Gy coordiante for Generator point
pub const G_Y_BYTES: [u8;32] = hex_literal::hex!("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8");

/// Represents a point P(x,y) on the elliptic curve
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EcPoint {
    Infinity,
    Point { x: FieldElement, y: FieldElement },
}

impl EcPoint {
    pub(crate) fn new(x: FieldElement, y: FieldElement) -> Self {
        EcPoint::Point { x, y }
    }

    pub(crate) fn add(self, other: Self) -> Self {
        match (self, other) {
            (EcPoint::Infinity, _) => other,
            (_, EcPoint::Infinity) => self,
            (EcPoint::Point { x: x1, y: y1 }, EcPoint::Point { x: x2, y: y2 }) => {
                // Case 1 if x1 == x2
                if x1 == x2 {
                    // Case 1a
                    // This means it has the same value of x but different values of y
                    // This means it is a vertical line and does not intersect at any point
                    if y1 != y2 {
                        return EcPoint::Infinity;
                    }
                    // Case 1b (y1 or y2 == 0)
                    // This means it is a vertical line and intersects at only one point
                    // Tangent is zero
                    if y1.value == U256::zero() || y2.value == U256::zero() {
                        return EcPoint::Infinity;
                    }
                    // Case 1c (y1 == y2)
                    // Point doubling
                    // If it has the same x and y for the 2 points
                    // P + P = 2P
                    // s(slope / differentiaton) = (3x^2 + a)/ 2y
                    let numerator = FieldElement::new(
                        multiply(U256::from(3), multiply(x1.value, x1.value)) + A,
                    );
                    let denominator = FieldElement::new(multiply(U256::from(2), y1.value));
                    //@note: This is where the division occurs, we try to avoid this here
                    let s = numerator / denominator;
                    let s_squared = FieldElement::new(multiply(s.value, s.value));
                    let x3 = s_squared - x1 - x2;
                    let y3 = FieldElement::new(multiply(s.value, (x1 - x3).value)) - y1;
                    return EcPoint::Point { x: x3, y: y3 };
                } else {
                    // Case 2 (x1 != x2)
                    // Point Addition (P + Q where P!=Q)
                    // s = (y2-y1)/(x2-x1)
                    //@note: This is where the division occurs, we try to avoid this here
                    let s = (y2 - y1) / (x2 - x1);
                    let s_squared = FieldElement::new(multiply(s.value, s.value));
                    let x3 = s_squared - x1 - x2;
                    let y3 = FieldElement::new(multiply(s.value, (x1 - x3).value)) - y1;
                    return EcPoint::Point { x: x3, y: y3 };
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::P;
    use super::*;

    /// Helper function to verify a point is on the curve: y^2 = x^3 + 7 (mod P)
    fn is_on_curve(point: &EcPoint) -> bool {
        match point {
            EcPoint::Infinity => true,
            EcPoint::Point { x, y } => {
                let y_squared = multiply(y.value, y.value);
                let x_cubed = multiply(multiply(x.value, x.value), x.value);
                let right_side = FieldElement::new(x_cubed + B);
                y_squared == right_side.value
            }
        }
    }

    /// Helper to get the secp256k1 generator point G
    fn get_generator() -> EcPoint {
        let gx = U256::from_big_endian(&G_X_BYTES);
        let gy = U256::from_big_endian(&G_Y_BYTES);

        EcPoint::Point {
            x: FieldElement::new(gx),
            y: FieldElement::new(gy),
        }
    }

    #[test]
    fn test_identity_property_infinity_left() {
        // Test: O + P = P (where O is the point at infinity)
        let g = get_generator();

        let result = EcPoint::Infinity.add(g);
        assert_eq!(result, g);
    }

    #[test]
    fn test_identity_property_infinity_right() {
        // Test: P + O = P (where O is the point at infinity)
        let g = get_generator();

        let result = g.add(EcPoint::Infinity);
        assert_eq!(result, g);
    }

    #[test]
    fn test_commutativity() {
        // Test: P + Q = Q + P
        // Using G and 2G as our test points (both are on the curve)
        let g = get_generator();
        let two_g = g.add(g);

        let p_plus_q = g.add(two_g);
        let q_plus_p = two_g.add(g);

        assert_eq!(p_plus_q, q_plus_p);
    }

    #[test]
    fn test_associativity() {
        // Test: (P + Q) + R = P + (Q + R)
        // Using G, 2G, and 3G as our test points
        let g = get_generator();
        let two_g = g.add(g);
        let three_g = g.add(two_g);

        let left = g.add(two_g).add(three_g);
        let right = g.add(two_g.add(three_g));

        assert_eq!(left, right);
    }

    #[test]
    fn test_point_doubling() {
        // Test: P + P = 2P (point doubling)
        let g = get_generator();

        let doubled = g.add(g);

        // Verify the result is not infinity
        assert_ne!(doubled, EcPoint::Infinity);

        // Verify the doubled point is on the curve
        assert!(is_on_curve(&doubled));
    }

    #[test]
    fn test_point_doubling_at_zero_y() {
        // Test: P + P = O when y = 0 (tangent is vertical)
        // We need to find a point where y = 0
        // For y^2 = x^3 + 7, when y = 0: 0 = x^3 + 7 (mod P)
        // This means x^3 = -7 (mod P) = P - 7

        // For this test, we'll create a point with y = 0
        // Note: This point may not actually be on the curve
        let x = FieldElement::new(U256::from(5));
        let y = FieldElement::new(U256::zero());
        let p = EcPoint::Point{ x, y };

        let result = p.add(p);
        assert_eq!(result, EcPoint::Infinity);
    }

    #[test]
    fn test_inverse_property() {
        // Test: P + (-P) = O (where -P has same x but negated y)
        let g = get_generator();

        // Extract the coordinates of G
        if let EcPoint::Point{ x, y } = g {
            // -G has the same x coordinate but y is negated (P - y in the field)
            let neg_y = FieldElement::new(P - y.value);
            let neg_g = EcPoint::Point { x, y: neg_y };

            let result = g.add(neg_g);
            assert_eq!(result, EcPoint::Infinity);
        } else {
            panic!("Generator should not be infinity");
        }
    }

    #[test]
    fn test_different_x_coordinates() {
        // Test point addition with different x coordinates
        // Using G and 2G which have different x coordinates
        let g = get_generator();
        let two_g = g.add(g);

        let result = g.add(two_g);

        // Result should not be infinity for different x coordinates
        assert_ne!(result, EcPoint::Infinity);

        // Verify the result is on the curve
        assert!(is_on_curve(&result));
    }

    #[test]
    fn test_secp256k1_generator_point_doubling() {
        // Test with a known point on secp256k1
        // Using the generator point G of secp256k1
        let g = get_generator();

        // Verify G is on the curve
        assert!(is_on_curve(&g));

        // Double the generator point: 2G
        let two_g = g.add(g);

        // Verify 2G is on the curve
        assert!(is_on_curve(&two_g));

        // Verify 2G is not infinity
        assert_ne!(two_g, EcPoint::Infinity);
    }

    #[test]
    fn test_secp256k1_generator_point_addition() {
        // Test G + 2G = 3G
        let g = get_generator();

        let two_g = g.add(g);
        let three_g = g.add(two_g);

        // Verify 3G is on the curve
        assert!(is_on_curve(&three_g));

        // Verify 3G is not infinity
        assert_ne!(three_g, EcPoint::Infinity);

        // Also verify that 2G + G = G + 2G (commutativity)
        let three_g_alt = two_g.add(g);
        assert_eq!(three_g, three_g_alt);
    }

    #[test]
    fn test_multiple_additions() {
        // Test: 2(P + Q) = 2P + 2Q
        // Using G and 2G
        let g = get_generator();
        let two_g = g.add(g);

        // Left side: 2(G + 2G) = 2(3G) = 6G
        let three_g = g.add(two_g);
        let six_g_left = three_g.add(three_g);

        // Right side: 2G + 4G = 6G
        let four_g = two_g.add(two_g);
        let six_g_right = two_g.add(four_g);

        assert_eq!(six_g_left, six_g_right);
    }

    #[test]
    fn test_infinity_plus_infinity() {
        // Test: O + O = O
        let result = EcPoint::Infinity.add(EcPoint::Infinity);
        assert_eq!(result, EcPoint::Infinity);
    }

    #[test]
    fn test_point_on_curve_validation() {
        // Test that our helper function correctly validates points on the curve

        // Point at infinity is always on the curve
        assert!(is_on_curve(&EcPoint::Infinity));

        // Test with the generator point (known to be on curve)
        let g = get_generator();
        assert!(is_on_curve(&g));
    }

    #[test]
    fn test_addition_preserves_curve_property() {
        // Test that adding two points on the curve results in a point on the curve
        let g = get_generator();

        // Start with G (on curve)
        assert!(is_on_curve(&g));

        // 2G should be on curve
        let two_g = g.add(g);
        assert!(is_on_curve(&two_g));

        // 3G should be on curve
        let three_g = g.add(two_g);
        assert!(is_on_curve(&three_g));

        // 4G should be on curve
        let four_g = two_g.add(two_g);
        assert!(is_on_curve(&four_g));

        // 5G should be on curve
        let five_g = two_g.add(three_g);
        assert!(is_on_curve(&five_g));
    }

    #[test]
    fn test_known_point_doubling_result() {
        // Test that 2G produces the expected result
        // This is a regression test with known values
        let g = get_generator();
        let two_g = g.add(g);

        // The expected coordinates of 2G on secp256k1
        let expected_x = U256::from_dec_str(
            "89565891926547004231252920425935692360644145829622209833684329913297188986597",
        )
        .unwrap();
        let expected_y = U256::from_dec_str(
            "12158399299693830322967808612713398636155367887041628176798871954788371653930",
        )
        .unwrap();

        let expected_two_g = EcPoint::Point {
            x: FieldElement::new(expected_x),
            y: FieldElement::new(expected_y),
        };

        assert_eq!(two_g, expected_two_g);
    }

    #[test]
    fn test_distributive_property() {
        // Test: k(P + Q) = kP + kQ for k = 2
        // This is: 2(G + 2G) = 2G + 4G
        let g = get_generator();
        let two_g = g.add(g);

        // Left side: 2(G + 2G) = 2(3G) = 6G
        let three_g = g.add(two_g);
        let left = three_g.add(three_g);

        // Right side: 2G + 4G = 6G
        let four_g = two_g.add(two_g);
        let right = two_g.add(four_g);

        assert_eq!(left, right);
    }

    #[test]
    fn test_repeated_doubling() {
        // Test: 2(2(2G)) = 8G
        let g = get_generator();

        let two_g = g.add(g);
        let four_g = two_g.add(two_g);
        let eight_g = four_g.add(four_g);

        // All should be on the curve
        assert!(is_on_curve(&two_g));
        assert!(is_on_curve(&four_g));
        assert!(is_on_curve(&eight_g));

        // None should be infinity
        assert_ne!(two_g, EcPoint::Infinity);
        assert_ne!(four_g, EcPoint::Infinity);
        assert_ne!(eight_g, EcPoint::Infinity);
    }
}
