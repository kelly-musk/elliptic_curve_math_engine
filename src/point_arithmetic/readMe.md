# Point arithmetic

## File [modular_arithmetic.rs](modular_arithmetic.rs)
- Implements field arithmetic operations for secp256k1
- Includes Add, Sub, Mul, Div operations for FieldElement

## File [point.rs](point.rs)
- Implements elliptic curve point addition
- Includes comprehensive tests for point arithmetic properties
### Features: 
- Data serialization (Public Keys, Signatures)
- 33 Bytes (32 bytes $x$ + 1 byte sign)
- Expensive (Requires 1 Inversion)
- Unique flag or specific logic

## File [projective_point.rs](projective_point.rs)
- Implements projective point addition
- Includes comprehensive tests for projective point arithmetic properties
### Features: 
- Runtime Computation (Verifying, Signing)
- 96 Bytes (Three 32-byte integers)
- Cheap (12 Multiplications, 0 Inversions)
- Contains Explicit $Z=0$ for inifinity flag


# Flow
1. Collect value in [EcPoint](point.rs)
2. Convert to [JacobianPoint](projective_point.rs)
3. Perform operations
4. Convert back to [EcPoint](point.rs)
