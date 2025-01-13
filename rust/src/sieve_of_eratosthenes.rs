pub fn get_primes_up_to(n: usize) -> Vec<usize> {
    if n < 2 {
        return vec![]; // No primes less than 2
    }

    // Create a boolean vector to track prime numbers
    let mut is_prime = vec![true; n + 1];
    is_prime[0] = false ; // 0 is not prime
    is_prime[1] = false; // 1 is not prime

    // Sieve of Eratosthenes
    let sqrt_n = (n as f64).sqrt() as usize;
    for i in 2..=sqrt_n {
        if is_prime[i] {
            for j in (i * i..=n).step_by(i) {
                is_prime[j] = false; // Mark multiples of i as non-prime
            }
        }
    }

    // Collect all prime numbers into a vector
    is_prime
        .iter()
        .enumerate()
        .filter(|&(_, &prime)| prime)
        .map(|(i, _)| i)
        .collect()
}
    