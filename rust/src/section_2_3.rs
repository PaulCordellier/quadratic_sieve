use bnum::{cast::As, types::U512};
// Possible other packages for long numbers:
// https://github.com/rust-num/num-bigint?tab=readme-ov-file#alternatives

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

#[allow(non_snake_case)]
pub fn algo_2_3_8(p: u128, a: u128) -> u32 {

    // In this function, I will be using u32 because the exp function requires u32 numbers.

    if jacobi_symbol(a, p) != 1 {
        panic!("The jacobi symbol is equal to {}, a = {}, p = {}", jacobi_symbol(a, p), a, p);
    }

    match p % 8 {
        3 | 7 => {
            let a: U512 = (a % p).as_();
            let exponent: u32 = (p as u32 + 1) / 4;
            let p: U512 = p.as_();
            
            return (a.pow(exponent) % p).as_::<u32>();
        },
        5 => {
            let a: u128 = a % p;
            let p_as_U512: U512 = p.as_();
            let mut x: U512 = a.as_::<U512>().pow((p as u32 + 3) / 8) % p_as_U512;
            let c: U512 = x.pow(2) % p_as_U512;

            if c != a.into() {
                x = (x * (2.as_::<U512>().pow((p as u32 - 1) / 4))) % p_as_U512;
            }

            return x.as_::<u32>();
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

            let a: u128 = a % p;
            let d: u128 = d % p;

            let A: u128 = a.pow(t) % p;
            let D: u128 = d.pow(t) % p;

            let mut m: u32 = 0;

            let A_as_U512: U512 = A.as_::<U512>();
            let D_as_U512: U512 = D.as_::<U512>();
            let p_as_U512: U512 = p.as_::<U512>();
            
            for i in 0..s {
                
                let exponent = 2_u32.pow(s - 1 - i);
                
                // I am deeply sorry to anyone reading this.
                if ((A_as_U512 * D_as_U512.pow(m) % p_as_U512).pow(exponent) % p_as_U512) % p_as_U512 == p_as_U512 - U512::ONE {
                    m += 2_u32.pow(i);
                }
            }

            return (a.pow((t + 1) / 2) * D.pow(m / 2) % p) as u32;
        }
        _ => panic!("The number p = {} should be even.", p)
    }
}