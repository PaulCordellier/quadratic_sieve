import math
import sieve_of_eratosthenes as eratos
import section_2_3
import get_zero_vector_combination as get_zero_vector


# Basic implementation of the quadratic sieve (page 276):
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


    ## 2. Sieving:

    S: list[tuple[int, int]] = []
    vector_matrix = []     # this is the matrix used in step 4
    vector_matrix_mod_2 = []     # this is the matrix used in step 3

    # As written in remark (2) of page 277, we won't sieve for small primes.
    primes_for_sieving = [prime for prime in primes if prime > 5]

    roots_a_for_sieving = roots_a[-len(primes_for_sieving):]

    prime_logarithms_for_sieving = [math.ceil(math.log2(x)) for x in primes_for_sieving]

    start_of_number_range = math.ceil(math.sqrt(n))

    number_range_size = 500
    error_margin = math.ceil(math.log2(B))

    number_of_trails = 0

    right_vector_matrix_is_found = False

    nb_of_tested_numbers = 0
    nb_of_sieved_numbers = 0
    nb_of_b_smooth_found = 0
    

    while not right_vector_matrix_is_found:

        number_range = list(range(start_of_number_range, start_of_number_range + number_range_size))

        number_of_trails += 1

        sum_of_logarithms = [0] * number_range_size

        # Sieving the number range here:

        for i in range(len(primes_for_sieving)):

            prime = primes_for_sieving[i]
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

            nb_of_tested_numbers += 1

            # if the sieve doesn't determine the number as b-smooth (this is approximate)
            if not sum_of_logarithms[i] >= math.log2(number_range[i]) - error_margin:
                continue

            nb_of_sieved_numbers += 1

            x = number_range[i]
            x2_n = x ** 2 - n     # x2_n = x² − n

            # So here the number x is probably B-smooth. Because the sieving with logarithm
            # isn't exact, we check here that the number is really B-smooth and make the
            # vector v(x² − n) of step 3 at the same time.

            # Establishing prime factorization of x2_n:
            returned_vectors = prime_factorization(x2_n, primes)

            if returned_vectors == None:    # In that case, the number isn't b-smooth
                continue

            nb_of_b_smooth_found += 1

            (vector, vector_mod_2) = returned_vectors

            # Here, the number is b-smooth for sure!
            S.append((x, x2_n))
            vector_matrix.append(vector)
            vector_matrix_mod_2.append(vector_mod_2)

            if len(S) < K + 10:
                continue

            vector_indexes = get_zero_vector.get_zero_vector_combination(vector_matrix_mod_2, K)

            if vector_indexes == None:
                continue

            # Here, if we add all of the vectors represented by the vector indexes mod 2,
            # the result is the zero vector (the zero vector needs to be obtained in step 3)

            S = [S[j] for j in vector_indexes]
            vector_matrix = [vector_matrix[j] for j in vector_indexes]
            vector_matrix_mod_2 = [vector_matrix_mod_2[j] for j in vector_indexes]

            right_vector_matrix_is_found = True
            break

        if number_of_trails > 1000:
            print("too much trails !")
            return 0
    
        start_of_number_range += number_range_size


    print("nb_of_tested_numbers : " + str(nb_of_tested_numbers))
    print("nb_of_sieved_numbers : " + str(nb_of_sieved_numbers))
    print("nb_of_b_smooth_found : " + str(nb_of_b_smooth_found))


    ## 3. Linear algebra:

    # Establishing the matrices has been done in step 2.
    # The addition of all vectors in vector_matrix_mod_2
    # will result in a zero vector.


    ## 4. Factorization:
    
    x = 1
    y = 1
    vector_for_y = [0] * K

    for i in range(len(vector_matrix)):

        x = (x + S[i][1]) % n

        vector_for_y = [vector_for_y[j] + vector_matrix[i][j] for j in range(K)]

    for i in range(K):
        y = y * pow(primes[i], vector_for_y[i] // 2, n) % n

    return math.gcd(x - y, n)


# Function for step 1:
def L(n: int):
    return math.e ** math.sqrt(math.log(n) * math.log(math.log(n)))


# Function for step 2:
def prime_factorization(n: int, primes: list[int]):
    """
    Returns an array for prime exponents and the same array mod 2,
    depending on the given primes
    """

    vector = [0] * len(primes)
    vector_mod_2 = [0] * len(primes)

    for i in range(len(primes)):

        while n / primes[i] % 1 == 0:

            n = n // primes[i]

            vector[i] += 1
            vector_mod_2[i] = 1 if vector_mod_2[i] == 0 else 0

    if n != 1:
        return None

    return (vector, vector_mod_2)



result_zero_or_one = 0
result_with_right_guess = 0
result_with_wrong_guess = 0

for i in range(157289988114331, 157289988114459, 2):
    result = quadratic_sieve(i)

    print(f"result for {i}: {result}")

    if result not in [0, 1]:
        if i % result == 0:
            print("Verification succeded !")
            result_with_right_guess += 1
        else:
            print("Verification failed :(")
            result_with_wrong_guess += 1
    else:
        result_zero_or_one += 1

    print()

print()
print("number of function results with :")
print(f"    zero or one : {result_zero_or_one}")
print(f"    right guess : {result_with_right_guess}")
print(f"    wrong guess : {result_with_wrong_guess}")