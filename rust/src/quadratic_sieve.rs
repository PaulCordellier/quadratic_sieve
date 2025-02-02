use std::usize;
use std::vec;
use num_integer::Integer;
use num_bigint::BigUint;

use crate::get_zero_vector_combination::get_zero_vector_combination;
use crate::section_2_3;
use crate::sieve_of_eratosthenes;

#[allow(non_snake_case)]
pub fn quadratic_sieve(n: u128, sieve_error_margin: u32) -> Result<QuadraticSieveAlogResults, &'static str> {

    // 1. Initialisation

    let B = L(n).isqrt();

    let primes: Vec<usize> = sieve_of_eratosthenes::get_primes_up_to(B)
                                .iter()
                                .filter(|&&prime| prime == 2 || section_2_3::jacobi_symbol(n, prime as u128) == 1)
                                .map(|prime| *prime as usize)
                                .collect();

    let K = primes.len();

    let roots_a: Vec<usize> = primes.iter()
                                    .map(|&prime| {
                                        if prime == 2 {
                                            1
                                        } else {
                                            section_2_3::algo_2_3_8(prime as u128, n) as usize
                                        }
                                    })
                                    .collect();



    // 2. Sieving

    let mut S: Vec<(u128, u128)> = vec![];
    let mut vector_matrix: Vec<Vec<u32>> = vec![];          // this is the matrix used in step 4
    let mut vector_matrix_mod_2: Vec<Vec<u32>> = vec![];    // this is the matrix used in step 3

    let mut primes_for_sieving: Vec<usize> = primes.clone();

    const MINIMUM_LOGARITHM: u32 = 4;

    for i in 0..K {

        if primes_for_sieving[i].ilog2() + 1 > MINIMUM_LOGARITHM {
            break;
        }

        while primes_for_sieving[i].ilog2() + 1 <= MINIMUM_LOGARITHM {
            primes_for_sieving[i] *= primes[i];
        }
    }


    let prime_logarithms_for_sieving: Vec<u32> = primes_for_sieving.iter()
                                                                   .map(|&prime| prime.ilog2() + 1)
                                                                   .collect();

    let mut start_of_number_range: usize = (n as f32).sqrt().ceil() as usize + 1;

    const NUMBER_RANGE_SIZE: usize = 1000;

    let mut number_of_trails: usize = 0;

    let mut quadratic_sieve_alog_results = QuadraticSieveAlogResults {
        non_trivial_factor: 0,
        nb_of_tested_numbers: 0,
        nb_of_sieved_numbers: 0,
        nb_of_b_smooth_found: 0,
        nb_of_trails_to_find_zero_vector: 0
    };
    
    'while_zero_vector_combination_is_not_found: loop {
    
        number_of_trails += 1;

        let mut sum_of_logarithms: [u32; NUMBER_RANGE_SIZE] = [0; NUMBER_RANGE_SIZE];

        // For each prime:
        for i in 0..primes_for_sieving.len() {

            let prime = primes_for_sieving[i];
            let root_a = roots_a[i];

            // first_prime_occurence represents here the position of the first number of 
            // the number range where number ≡ 0 (mod prime)
            let mut first_prime_occurence = prime - start_of_number_range % prime;

            if first_prime_occurence == prime {
                first_prime_occurence = 0;
            }

            let prime_multiplicator_for_root_a = root_a / prime;

            // first_occurence_for_positive_a represents here the position of the first number of 
            // the number range where number ≡ root_a (mod prime) 
            let first_occurence_for_positive_a = first_prime_occurence + prime * (prime_multiplicator_for_root_a + 1) - root_a;

            // first_occurence_for_negative_a is the same as first_occurence_for_positive_a, but the root_a
            // is negative. The equation becomes number ≡ -root_a (mod prime)
            let first_occurence_for_negative_a = first_prime_occurence + root_a - prime * prime_multiplicator_for_root_a;

            for j in (first_occurence_for_positive_a..NUMBER_RANGE_SIZE).step_by(prime) {
                sum_of_logarithms[j] += prime_logarithms_for_sieving[i]
            }

            for j in (first_occurence_for_negative_a..NUMBER_RANGE_SIZE).step_by(prime) {
                sum_of_logarithms[j] += prime_logarithms_for_sieving[i]
            }
        }

        let mut threshold_to_recognize_b_smooth = ((start_of_number_range as u128).pow(2) - n).ilog2();
        threshold_to_recognize_b_smooth -= sieve_error_margin;
    
        // For each number in the number range:
        for i in 0..NUMBER_RANGE_SIZE {

            if sum_of_logarithms[i] < threshold_to_recognize_b_smooth {
                continue;
            }
            
            quadratic_sieve_alog_results.nb_of_sieved_numbers += 1;

            let x = (start_of_number_range + i) as u128;

            // Here, x has been sieved and is probably B-smooth.
            // We need to be sure though. The function prime_factorization returns None if
            // x isn't B-smooth. If it is, the vector of exponents of primes factors of x
            // is returned in vector. (The same vector (mod 2) is put in vector_mod_2).

            let x2_n = x.pow(2) - n;     // x2_n = x² − n

            let vector: Vec<u32>;
            let vector_mod_2: Vec<u32>;

            match prime_factorization(x2_n, &primes) {
                Some(x) => {
                    (vector, vector_mod_2) = x;
                },
                None => continue
            }

            quadratic_sieve_alog_results.nb_of_b_smooth_found += 1;

            S.push((x, x2_n));
            vector_matrix.push(vector);
            vector_matrix_mod_2.push(vector_mod_2);

            if S.len() < K + 10 {
                continue;
            }

            let vector_indexes: Vec<usize>;

            quadratic_sieve_alog_results.nb_of_trails_to_find_zero_vector += 1;

            // If a zero vector combination is found, we put the indexes of the vectors
            // that need to be added to get the  in the variable vector_indexes.
            // If no combination is found, we just continue the loop.
            match get_zero_vector_combination(&vector_matrix_mod_2, K) {
                Some(x) => vector_indexes = x,
                None => continue
            };

            // Here, if we add all of the vectors represented by the vector indexes mod 2,
            // the result is the zero vector (the zero vector needs to be obtained in step 3)
            
            assert_eq!(S.len(), vector_matrix.len());
            assert_eq!(S.len(), vector_matrix_mod_2.len());

            for i in (0..S.len()).rev() {
                if !vector_indexes.contains(&i) {
                    S.swap_remove(i);
                    vector_matrix.swap_remove(i);
                    vector_matrix_mod_2.swap_remove(i);
                }
            }

            quadratic_sieve_alog_results.nb_of_tested_numbers = (number_of_trails * NUMBER_RANGE_SIZE) as u64;

            break 'while_zero_vector_combination_is_not_found;
            // This break statement exits the two loops we are currently in
        }

        start_of_number_range += NUMBER_RANGE_SIZE;

        if number_of_trails >= 5000 {
            return Err("The algorithm didn't find enough B-smooth values");
        }
    }


    // 3. Linear algebra:

    // Establishing the matrices has been done in step 2.
    // The addition of all vectors in vector_matrix_mod_2
    // will result in a zero vector.


    // 4. Factorization:

    let mut x: u128 = 1;
    let mut vector_for_y = vec![0; K];
    
    for (i, vector) in vector_matrix.iter().enumerate() {

        x = (x + S[i].0) % n;

        for j in 0..K {
            vector_for_y[j] += vector[j];
        }

    }

    let mut y: BigUint = BigUint::from(1_u8);
    let n_as_biguint: BigUint = BigUint::from(n);
    
    for i in 0..K {
        let prime: BigUint = BigUint::from(primes[i]);
        let prime_exponent: BigUint = BigUint::from(vector_for_y[i] / 2);

        y *= prime.modpow(&prime_exponent, &n_as_biguint);
        y %= &n_as_biguint;
    }
    
    let y: u128 = y.try_into().expect("Couldn't convert the biguint y to an u128.");

    quadratic_sieve_alog_results.non_trivial_factor = n.gcd(&(y - x));

    // algos to use gcd :
    // - https://docs.rs/num-integer/0.1.41/num_integer/fn.gcd.html
    // - https://docs.rs/bnum/latest/bnum/?search=gcd
    
    return Ok(quadratic_sieve_alog_results);
    
}

pub struct QuadraticSieveAlogResults {
    pub non_trivial_factor: u128,
    pub nb_of_tested_numbers: u64,
    pub nb_of_sieved_numbers: u64,
    pub nb_of_b_smooth_found: u64,
    pub nb_of_trails_to_find_zero_vector: u64
}


/// Function used in step 3 of the basic quadratic sieve.
/// 
/// Returns None if the given primes can't factor n.
/// If the primes can factor n, the first given vector is a vector exponents of
/// the given primes factors.
/// The second given vector is the same than the first but (mod 2).
fn prime_factorization(n: u128, primes: &Vec<usize>) -> Option<(Vec<u32>, Vec<u32>)> {

    let mut vector = vec![0; primes.len()];
    let mut vector_mod_2 = vec![0; primes.len()];

    let mut n = n;

    for (i, &prime) in primes.iter().enumerate() {
        
        let prime = prime as u128;

        while n % prime == 0 {
            n /= prime;
            vector[i] += 1;
            vector_mod_2[i] = if vector_mod_2[i] == 0 { 1 } else { 0 };
        }
    }

    if n != 1 {
        return None;
    }

    Some((vector, vector_mod_2))

}

// Function for step 1:
#[allow(non_snake_case)]
fn L(n: u128) -> usize {
    let n = n as f64;
    (std::f64::consts::E).powf((n.ln() * n.ln().ln()).sqrt()) as usize
}