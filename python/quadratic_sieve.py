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

    print("odd_primes: " + str(primes))
    print("roots_a: " + str(roots_a))


    ## 2. Sieving:

    S = []
    matrix = []     # this is the matrix used in step 3

    # As written in remark (2) of page 277, we won't sieve for small primes.
    primes_for_sieving = [prime for prime in primes if prime > 5]

    roots_a_for_sieving = roots_a[-len(primes_for_sieving):]

    prime_logarithms_for_sieving = [int(math.log2(x)) for x in primes_for_sieving]

    start_of_next_number_range = math.ceil(math.sqrt(n))

    number_range_size = 500
    error_margin = math.ceil(math.log2(B))

    while len(S) < K + 1:

        number_range = list(range(start_of_next_number_range, start_of_next_number_range + number_range_size))
        start_of_next_number_range += number_range_size

        sum_of_logarithms = [0] * number_range_size

        # Sieving the number range here:

        for i in range(len(primes_for_sieving)):

            prime = primes_for_sieving[i]

            start_of_number_range = number_range[0]
            root_a = roots_a_for_sieving[i]

            # first_occurence here represents the position of the first number of number_range for which
            # number ≡ 0 (mod prime) 
            first_occurence = math.ceil(start_of_number_range / prime) * prime - start_of_number_range
            
            # first_occurence_for_positive_a represents here the position of the first number of 
            # number_range for which number ≡ root_a (mod prime) 
            if first_occurence + root_a >= prime:
                first_occurence_for_positive_a = first_occurence - prime + root_a
            else:
                first_occurence_for_positive_a = first_occurence + root_a

            for j in range(first_occurence_for_positive_a, number_range_size, prime):
                sum_of_logarithms[j] += prime_logarithms_for_sieving[i]

            # first_occurence_for_negative_a is the same as first_occurence_for_positive_a, but the root_a
            # is negative. The equation becomes number ≡ -root_a (mod prime) 
            if first_occurence - root_a < 0:
                first_occurence_for_negative_a = first_occurence + prime - root_a
            else:
                first_occurence_for_negative_a = first_occurence - root_a

            for j in range(first_occurence_for_negative_a, number_range_size, prime):
                sum_of_logarithms[j] += prime_logarithms_for_sieving[i]
                
        for i in range(number_range_size):
            if sum_of_logarithms[i] >= math.log2(number_range[i]) - error_margin:

                x = number_range[i]
                x2_n = x ** 2 - n     # x2_n = x² − n

                pair_to_add_in_S = (x, x2_n)

                # So here the number x is probably B-smooth. Because the sieving with logarithm
                # isn't exact, we check here that the number is really B-smooth and make the
                # vector v(x2 − n) of step 3 at the same time.

                # Establishing prime factorization of x2_n:

                vector = [0] * len(primes)

                for j in range(len(primes)):

                    while x2_n / primes[j] % 1 == 0:

                        x2_n = x2_n // primes[j]

                        vector[j] = 1 if vector[j] == 0 else 0
                
                if x2_n == 1:
                    # Here, the number is b-smooth for sure!
                    S.append(pair_to_add_in_S)
                    matrix.append(vector)

    S = S[:K+1]
    matrix = matrix[:K+1]

    print(S)
    
    
    ## 3. Linear algebra:

    # Establishing the matrix has been done in step 2.

    print(matrix)


    ## 4. Factorization:

    # ...


quadratic_sieve(158189981)