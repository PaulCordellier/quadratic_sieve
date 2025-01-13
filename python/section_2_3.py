# This code contains the functions we need from the section 2.3 of the book

def jacobi_symbol(a: int, m: int) -> int:
    """
    This is the algorithm 2.3.5, page 112 in the manual.
    Given positive odd integer m, and integer a, this algorithm
    returns the Jacobi symbol, which for m an odd prime is also
    the Legendre symbol.
    """

    if m % 2 == 0:
        raise ValueError("The integer m must be odd. m = " + m)
    
    if m < 0:
        raise ValueError("The integer m must be positive. m = "+ m)
    
    a = a % m

    t = 1

    while a != 0:
        while a % 2 == 0: # while a is even
            a = a / 2
            if m % 8 == 3 or m % 8 == 5:
                t = -t
        a, m = m, a
        if (a % 4 == 3 and m % 4 == 3):
            t = -t
        a = a % m

    if m == 1:
        return t
    
    return 0


def algo_2_3_8(p: int, a: int) -> int:
    """
    This algorithm is at the pdf page 114 in the manual.
    Given an odd prime p and an integer a with jacobi_symbol(a, p) = 1, 
    this algorithm returns a solution x to x^2 ≡ a (mod p).
    """

    if jacobi_symbol(a, p) != 1:
        raise ValueError(f"The jacobi symbol is equal to {jacobi_symbol(a, p)}, a = {a}, p = {p}")

    match p % 8:

        case 3 | 7:
            a = a % p
            return a ** ((p + 1) // 4) % p
        
        case 5:
            a = a % p
            
            x = a ** ((p + 3) // 8) % p

            c = (x ** 2) % p

            if c != a:
                x = (x * 2 ** ((p - 1) // 4)) % p

            return x
        
        case 1:

            ## Find a random integer d ∈ [2, p − 1] with jacobi_symbol(d, p) = -1:
            d = 0

            for i in range(2, p):
                if jacobi_symbol(i, p) == -1:
                    d = i
                    break

            if d == 0:
                raise ValueError("A random integer d wasn't found")

            ## Represent p − 1 = (2 ** s) * t, with t odd:

            s = 0
            t = p - 1

            while True:

                tmp = t / 2

                if tmp % 1 != 0:
                    # when the code comes to this point, t is odd   
                    break
    
                s += 1
                t = int(tmp)

            if p - 1 != (2 ** s) * t:
                raise ValueError("wrong s and t values")
            
            if t % 2 == 0:
                raise ValueError(f"the value t shouldn't be odd. t = {t}")
            
            A = pow(a, t, p)
            D = pow(d, t, p)

            m = 0

            for i in range(0, s):
                if (A * D ** m) ** (2 ** (s - 1 - i)) % p == p - 1:
                    m += 2 ** i

            return pow(a, (t + 1) // 2, p) * pow(D, m // 2, p) % p
        
        case _:
            raise ValueError("What? this number is even")
