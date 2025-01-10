import math
import sieve_of_eratosthenes as eratos

def get_b_smooth_numbers_with_multiplication(b: int, n: int) -> list[int]:
    """
    Returns every b-smooth numbers in the range [2, n]
    """

    if b < 2:
        raise ValueError("The condition 2 ≤ B should be true.")

    if b > math.sqrt(n):
        raise ValueError("The condition B ≤ √N should be true.")

    # Create the number range [2, N]:
    number_range = list(range(2, n + 1))

    # B-Smooth values of the numbers in the number range
    b_smooth_value_of_numbers = [1] * len(number_range)

    # Get all the prime numbers on [2, B]
    primes_numbers = eratos.get_primes_up_to(b)

    # For each prime number
    for prime_number in primes_numbers:
        prime_number_to_power = prime_number

        while prime_number_to_power <= n:
            # For each multiple of the prime_number_to_power
            for i in range(prime_number_to_power, n + 1, prime_number_to_power):
                # At this point, number_range[i - 2] is i and is a multiple of prime_number_to_power
                b_smooth_value_of_numbers[i - 2] *= prime_number

            prime_number_to_power *= prime_number

    b_smooth_numbers = []
    
    for i in range(len(number_range)):
        if b_smooth_value_of_numbers[i] == number_range[i]:
            b_smooth_numbers.append(number_range[i])

    return b_smooth_numbers


def get_b_smooth_numbers_with_logarithm(b: int, n: int) -> list[int]:
    """
    Returns approximately every b-smooth numbers in the range [2, n]
    """

    if b < 2:
        raise ValueError("The condition 2 ≤ B should be true.")

    if b > math.sqrt(n):
        raise ValueError("The condition B ≤ √N should be true.")

    # Create the number range [2, N]:
    number_range = list(range(2, n + 1))

    # Logarithmic sum values
    addition_of_logarithm = [0] * len(number_range)

    # Get all the prime numbers on [2, B]
    primes_numbers = eratos.get_primes_up_to(b)

    # The logarithm of the primes numbers
    prime_number_logarithm = [int(math.log2(x) + 0.5) for x in primes_numbers]

    # For each prime number
    for i in range(len(primes_numbers)):
        prime_number_to_power = primes_numbers[i]

        while prime_number_to_power <= n:
            # For each multiple of the prime_number_to_power
            for j in range(prime_number_to_power, n + 1, prime_number_to_power):
                # At this point, number_range[j - 2] is j and is a multiple of prime_number_to_power
                addition_of_logarithm[j - 2] += prime_number_logarithm[i]

            prime_number_to_power *= primes_numbers[i]

    error_margin = math.log2(b)

    b_smooth_numbers = []
    
    for i in range(len(number_range)):
        if addition_of_logarithm[i] >= math.log2(number_range[i]) - error_margin:
            b_smooth_numbers.append(number_range[i])

    return b_smooth_numbers
