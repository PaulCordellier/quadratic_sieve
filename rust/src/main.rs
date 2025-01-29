use std::env;

mod quadratic_sieve;
mod sieve_of_eratosthenes;
mod section_2_3;
mod get_zero_vector_combination;

fn main() {

    env::set_var("RUST_BACKTRACE", "1");

    let result = quadratic_sieve::quadratic_sieve(157289988114349);

    match result {
        Ok(algo_results) => {
            println!("nb_of_tested_numbers = {}", algo_results.nb_of_tested_numbers);
            println!("nb_of_sieved_numbers = {}", algo_results.nb_of_sieved_numbers);
            println!("nb_of_b_smooth_found = {}", algo_results.nb_of_b_smooth_found);
        },
        Err(error_message) => println!("{}", error_message) 
    }
}

