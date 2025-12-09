//! Deriving a keypair from points on the sec256pk1 curve
//! 
//! With the wierstrass formula
use super::{
    private_key::PrivateKey,
    pubkey::PublicKey
};
use primitive_types::U256;

use crate::point_arithmetic::{
    P, get_generator_affine,get_generator_jacobian
};
pub(crate) const CURVE_ORDER_HEX: &str = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";
pub struct KeyPair {
    /// scalar k
    pub private_key: PrivateKey,
    /// EcPoint (x,y)
    pub public_key: PublicKey,
}

impl KeyPair{
    pub fn generate() -> Self {
        todo!()
    }
}