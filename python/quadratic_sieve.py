import math
import sieve_of_eratosthenes as eratos
import section_2_3

# Functions:

def L(n: int):
    return math.e ** math.sqrt(math.log(n) * math.log(math.log(n)))

def quadratic_sieve(n: int):
    """
    Implements the quadratic sieve. It's basic implementation is
    at the pdf page 276 of the book.
    """

    ## 1. Initialization:

    B = math.ceil(L(n) ** 0.5)

    primes = eratos.get_primes_up_to(B)[1:]     # [1:] is used so we don't use the even prime 2
    primes = [2] + list(filter(lambda prime: section_2_3.jacobi_symbol(n, prime) == 1, primes))

    K = len(primes)

    roots_a = [1]

    for i in range(1, K):
        roots_a.append(section_2_3.algo_2_3_8(primes[i], n))

    print("odd_primes " + str(primes))
    print("roots_a " + str(roots_a))


    ## 2. Sieving:

    prime_logarithms = [int(math.log2(x)) for x in primes]
    start_of_next_number_range = math.ceil(math.sqrt(n))

    S = []

    number_range_size = 500
    error_margin = math.ceil(math.log2(B) / 2)

    while len(S) <= K + 1:
        number_range = list(range(start_of_next_number_range, start_of_next_number_range + number_range_size))
        start_of_next_number_range += number_range_size

        addition_of_logarithms = [0] * number_range_size
        # Sieving un the number range here

        for i in range(len(primes)):

            for j in range(math.ceil(number_range[0] / primes[i]) * primes[i] - number_range[0], number_range_size, primes[i]):
            
                addition_of_logarithms[j] += prime_logarithms[i] 


        for i in range(len(number_range)):
            if addition_of_logarithms[i] >= math.log2(number_range[i]) - error_margin:

                S.append((number_range[i], number_range[i] ** 2 - n))

    print(S)

    ## 3. Linear algebra:

    # ...

    ## 4. Factorization:

    # ...


quadratic_sieve(19484388)