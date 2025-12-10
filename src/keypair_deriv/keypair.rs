//! Deriving a keypair from points on the sec256pk1 curve
//! 
//! With the wierstrass formula
use super::{
    private_key::PrivateKey,
    pubkey::PublicKey
};
use primitive_types::U256;

use crate::point_arithmetic::{
    JacobianPoint, P, get_generator_affine, get_generator_jacobian
};

/// How many points exist on the curve
/// 
/// CURVE_ORDER_HEX
/// 
/// 0 < k < N < P
pub(crate) const N: &str = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";

#[derive(Debug)]
pub struct KeyPair {
    /// scalar k
    pub private_key: PrivateKey,
    /// EcPoint (x,y)
    pub public_key: PublicKey,
}

impl KeyPair{
    pub fn generate() -> Self {
        // set the n to the N
        let n = U256::from_str_radix(N, 16).unwrap();
        // We will be accepting Affine / Ecpoint co-ordinates
        let g_affine = get_generator_affine();
        let g_jacobian = JacobianPoint::from(g_affine);
        todo!()
    }
}