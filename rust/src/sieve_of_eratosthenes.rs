pub fn get_primes_up_to(n: usize) -> Vec<u32> {

    if n < 2 {
        panic!("In the function get_primes_up_to(n) n shouldn't be less than 2. n = {}", n);
    }

    // Create a boolean vector to track prime numbers
    let mut is_prime: Vec<bool>  = vec![true; n + 1];
    is_prime[0] = false; // 0 is not prime
    is_prime[1] = false; // 1 is not prime

    // Sieve of Eratosthenes
    for i in 2..=n.isqrt() {
        if is_prime[i] {
            for j in (i * i..=n).step_by(i) {
                is_prime[j] = false; // Mark multiples of i as non-prime
            }
        }
    }

    is_prime.iter()
            .enumerate()
            .filter(|(_, is_prime)| **is_prime )
            .map(|(index, _)| index as u32)
            .collect()
}
    