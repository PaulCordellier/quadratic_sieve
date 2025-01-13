import math
import sieve_of_eratosthenes as eratos
import section_2_3

# Function for step 1:
def L(n: int):
    return math.e ** math.sqrt(math.log(n) * math.log(math.log(n)))

find_zero_vector_tests = 0

# Function for step 3:
def find_zero_vector(matrix: list[list[int]], vector_combination: list[bool], vector_addition: list[int] | None = None) -> tuple[list[bool], list[int]] | None:

    global find_zero_vector_tests
    find_zero_vector_tests += 1
    
    similatities_between_vectors = [0] * len(matrix)

    if vector_addition != None:
    
        for i in range(len(matrix)):

            if vector_combination[i] == True:
                continue
                
            similatities = sum(1 for x, y in zip(matrix[i], vector_addition) if x == y)

            similatities_between_vectors[i] = similatities


    best_indexes = sorted(range(len(similatities_between_vectors)), key=lambda i: similatities_between_vectors[i], reverse=True)

    print(vector_addition)
    # print()
    # print(similatities_between_vectors)
    # print(best_indexes)

    for i in best_indexes:
        
        if vector_combination[i] == True:
            continue

        check_this_vector = False

        if vector_addition == None:
            check_this_vector = True
        else:
            for j in range(len(matrix[i])):
                if matrix[i][j] == vector_addition[j] == 1:
                    check_this_vector = True
                    break

        if not check_this_vector:
            continue

        new_vector_combination = [vector_combination[j] if j != i else True for j in range(len(vector_combination))]

        if vector_addition == None:
            new_vector_addition = matrix[i]
        else:
            new_vector_addition = [0 if vector_addition[j] == matrix[i][j] else 1 for j in range(len(vector_addition))]


        # conunt the number of zeroes
        number_of_ones = sum(value for value in new_vector_addition)

        if number_of_ones == 0:
            return (new_vector_combination, new_vector_addition)
        
        if number_of_ones > len(new_vector_addition) * 0.2 + 2:
            return None

        return_values = find_zero_vector(matrix, new_vector_combination, new_vector_addition)

        if return_values == None:
            continue

        if all(value == 0 for value in return_values[1]):
            return return_values
            
    return None


# Basic implementation of the quadratic sieve (page 276):
def quadratic_sieve(n: int):
    """
    Implements the quadratic sieve. It's basic implementation is
    at the pdf page 276 of the book.
    """

    ## 1. Initialization:

    print("step 1")

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

    print("step 2")

    S: list[tuple[int, int]] = []
    vector_matrix = []     # this is the matrix used in step 4
    vector_matrix_mod_2 = []     # this is the matrix used in step 3

    # As written in remark (2) of page 277, we won't sieve for small primes.
    primes_for_sieving = [prime for prime in primes if prime > 5]

    roots_a_for_sieving = roots_a[-len(primes_for_sieving):]

    prime_logarithms_for_sieving = [int(math.log2(x)) for x in primes_for_sieving]

    start_of_next_number_range = math.ceil(math.sqrt(n))

    number_range_size = 500
    error_margin = math.ceil(math.log2(B))

    number_of_trails = 0

    while len(S) < K + 1:

        number_range = list(range(start_of_next_number_range, start_of_next_number_range + number_range_size))
        start_of_next_number_range += number_range_size

        number_of_trails += 1

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
                # vector v(x² − n) of step 3 at the same time.

                # Establishing prime factorization of x2_n:

                vector = [0] * len(primes)
                vector_mod_2 = [0] * len(primes)

                for j in range(len(primes)):

                    while x2_n / primes[j] % 1 == 0:

                        x2_n = x2_n // primes[j]

                        vector[j] += 1
                        vector_mod_2[j] = 1 if vector_mod_2[j] == 0 else 0
                
                if x2_n == 1:
                    # Here, the number is b-smooth for sure!
                    S.append(pair_to_add_in_S)
                    vector_matrix.append(vector)
                    vector_matrix_mod_2.append(vector_mod_2)

        if number_of_trails > 10000:
            return 0

    S = S[:K + 1]
    vector_matrix = vector_matrix[:K + 1]
    vector_matrix_mod_2 = vector_matrix_mod_2[:K + 1]

    print(S)
    
    
    ## 3. Linear algebra:

    print("step 3")

    # Establishing the matrix has been done in step 2.

    print(vector_matrix)
    print(vector_matrix_mod_2)

    print()

    returned_value = find_zero_vector(vector_matrix_mod_2, [False] * (K + 1))

    if returned_value == None:
        return 0

    vector_combination = returned_value[0]

    for i in range(len(vector_matrix_mod_2)):
        if not vector_combination[i]:
            continue
        print(vector_matrix_mod_2[i])


    ## 4. Factorization:

    print("step 4")

    x = 1
    y = 1
    vector_for_y = [0] * K

    print()

    for i in range(K + 1):

        if not vector_combination[i]:
            continue

        x = (x + S[i][0]) % n

        vector_for_y = [vector_for_y[j] + vector_matrix[i][j] for j in range(K)]

    print(vector_for_y)

    print()

    x = x % n

    for i in range(K):
        y = (y * (primes[i] ** (vector_for_y[i] // 2))) % n

    y = y

    print(f"x = {x}, y = {y}")
    print(f"find_zero_vector_tests = {find_zero_vector_tests}")

    return math.gcd(x - y, n)

# print(quadratic_sieve(1581899))
print(quadratic_sieve(158189981))
# print(quadratic_sieve(157289988331))
# print(quadratic_sieve(157289988114331))