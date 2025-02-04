use num_bigint::BigUint;
// Possible other packages for long numbers:
// https://github.com/rust-num/num-bigint?tab=readme-ov-file#alternatives

/// This is the algorithm 2.3.5, page 112 in the manual.
/// Given positive odd integer m, and integer a, this algorithm
/// returns the Jacobi symbol, which for m an odd prime is also
/// the Legendre symbol.
pub fn jacobi_symbol(mut a: u128, mut m: u128) -> i8 {

    if m % 2 == 0 {
        panic!("The integer m must be odd. m = {}", m);
    }

    a = a % m;
    let mut t: i8 = 1;

    while a != 0 {
        while a % 2 == 0 {  // while a is even
            a /= 2;
            if m % 8 == 3 || m % 8 == 5 {
                t = -t;
            }
        }
        std::mem::swap(&mut a, &mut m);
        if a % 4 == 3 && m % 4 == 3 {
            t = -t;
        }
        a = a % m;
    }

    if m == 1 {
        return t;
    }

    0
}

/// This algorithm is at the pdf page 114 in the manual.
/// Given an odd prime p and an integer a with jacobi_symbol(a, p) = 1, 
/// this algorithm returns a solution x to x^2 ≡ a (mod p).
#[allow(non_snake_case)]
pub fn algo_2_3_8(p: u128, a: u128) -> u32 {

    // In this function, I will be using u32 because the exp function requires u32 numbers.

    if jacobi_symbol(a, p) != 1 {
        panic!("The jacobi symbol is equal to {}, a = {}, p = {}", jacobi_symbol(a, p), a, p);
    }

    let biguint_1: &BigUint = &BigUint::from(1_u8);
    let biguint_2: &BigUint = &BigUint::from(2_u8);
    let biguint_3: &BigUint = &BigUint::from(3_u8);
    let biguint_4: &BigUint = &BigUint::from(4_u8);
    let biguint_8: &BigUint = &BigUint::from(8_u8);

    let result: BigUint = match p % 8 {
        3 | 7 => {
            let a = BigUint::from(a % p);
            let p = BigUint::from(p);
            let exponent = (&p + biguint_1) / biguint_4;

            a.modpow(&exponent, &p)
        },
        5 => {
            let a = BigUint::from(a % p);
            let p = BigUint::from(p);
            let mut x = a.modpow(&((&p + biguint_3) / biguint_8), &p);
            let c = x.modpow(biguint_2, &p);

            if c != a.into() {
                x = x * (biguint_2.modpow(&((&p - biguint_1) / biguint_4), &p));
            }

            x % p
        }
        1 => {
            // Find a random integer d ∈ [2, p − 1] with jacobi_symbol(d, p) = -1:
            let mut d = 0;
            for i in 2..p {
                if jacobi_symbol(i, p) == -1 {
                    d = i;
                    break;
                }
            }

            if d == 0 {
                panic!("A random integer d wasn't found");
            }

            // Represent p − 1 = (2 ** s) * t, with t odd:
            let mut s: u32 = 0;
            let mut t: u32 = p as u32 - 1;

            while t % 2 == 0 {
                t /= 2;
                s += 1;
            }

            if p as u32 - 1 != 2_u32.pow(s) * t || t % 2 == 0 {
                panic!("Incorrect values for s and/or t");
            }

            let p = BigUint::from(p);
            let t = BigUint::from(t);

            let a = BigUint::from(a) % &p;
            let d = BigUint::from(d) % &p;

            let A = a.modpow(&t, &p);
            let D = d.modpow(&t, &p);

            let mut m = BigUint::from(0_u8);

            for i in 0..s {
                
                let exponent = biguint_2.pow(s - 1 - i);
                
                if ((&A * D.modpow(&m, &p)).modpow(&exponent, &p)) % &p == &p - biguint_1 {
                    m += 2_u32.pow(i);
                }
            }

            a.modpow(&((t + biguint_1) / biguint_2), &p) * D.modpow(&(m / biguint_2), &p)
        }
        _ => panic!("The number p = {} should be even.", p)
    };

    biguint_to_u32(&result)
}

fn biguint_to_u32(biguint: &BigUint) -> u32 {
    let u32_digits = biguint.to_u32_digits();

    if u32_digits.len() == 0 {
        return 0;
    }

    biguint.to_u32_digits()[0]
}