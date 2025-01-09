def get_primes_up_to(n : int) -> list[int]:
    primes = []

    if n < 2:
        return primes  # No primes less than 2

    # Create a boolean list to track prime numbers
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False  # 0 and 1 are not prime

    # Sieve of Eratosthenes
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False  # Mark multiples of i as non-prime

    # Collect all prime numbers into a list
    for i in range(2, n + 1):
        if is_prime[i]:
            primes.append(i)

    return primes