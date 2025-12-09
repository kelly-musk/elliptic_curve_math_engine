//! Deriving a keypair from points on the sec256pk1 curve
//! 
//! With the wierstrass formula
use super::{
    private_key::PrivateKey,
    pubkey::PublicKey
};
pub struct KeyPair {
    pub private_key: PrivateKey,
    pub public_key: PublicKey,
}