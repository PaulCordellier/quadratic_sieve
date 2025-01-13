use std::io;

mod sieve_of_eratosthenes;

fn main() {

    // Création d'un buffer pour lire l'entrée utilisateur
    let mut input = String::new();

    // Lecture de l'entrée utilisateur
    io::stdin()
        .read_line(&mut input)
        .expect("IO error.");

    // Tentative de conversion de l'entrée en usize
    let n : usize = input.trim().parse::<usize>().expect("Come on write a positive number!");


    sieve_of_eratosthenes::get_primes_up_to(n);
    println!("Hello, world!");

    let array: [i8; 9] = [0; 9];

    let mut array2 = array.clone();

    array2[1] = 1;
}

