

use primitive_types::U256;

/// Private key must be a scalar k that satisfies 0 < k < n
/// 
/// n being the order / number of elements in the curve
#[derive(Debug)]
pub struct PrivateKey(pub U256);