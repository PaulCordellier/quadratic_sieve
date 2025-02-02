use std::time::{Duration, SystemTime};

mod quadratic_sieve;
mod sieve_of_eratosthenes;
mod section_2_3;
mod get_zero_vector_combination;

fn main() {

    const TEST_START_NUMBER: u128 = 10_u128.pow(20);
    const NUMBER_OF_TESTS: u128 = 100;
    const SIEVE_ERROR_MARGIN: u32 = 15;

    let mut number_of_right_guess = 0;
    let mut number_of_wrong_guess = 0;
    let mut number_of_error_cant_find_enough_b_smooth = 0;
    let mut number_of_tested_numbers = 0;
    let mut number_of_sieved_numbers = 0;
    let mut number_of_actual_b_smooth_numbers = 0;
    let mut number_of_trails_to_find_zero_vector = 0;

    let mut average_time_of_function = Duration::new(0, 0);

    for n in (TEST_START_NUMBER .. TEST_START_NUMBER + NUMBER_OF_TESTS * 2).step_by(2) {

        let time_before_quadratic_sieve = SystemTime::now();

        let result = quadratic_sieve::quadratic_sieve(n, SIEVE_ERROR_MARGIN);

        let elapsed_time = time_before_quadratic_sieve.elapsed().expect("System time error");

        average_time_of_function += elapsed_time;

        let algo_results = match result {
            Ok(algo_results) => algo_results,
            Err(_error_message) => {
                number_of_error_cant_find_enough_b_smooth += 1;
                continue;
            }
        };

        if algo_results.non_trivial_factor != 1 && n % algo_results.non_trivial_factor == 0 {
            number_of_right_guess += 1;
        } else {
            number_of_wrong_guess += 1;
            continue;
        }

        number_of_tested_numbers += algo_results.nb_of_tested_numbers;
        number_of_sieved_numbers += algo_results.nb_of_sieved_numbers;
        number_of_actual_b_smooth_numbers += algo_results.nb_of_b_smooth_found;
        number_of_trails_to_find_zero_vector += algo_results.nb_of_trails_to_find_zero_vector;
    }

    average_time_of_function /= NUMBER_OF_TESTS as u32;

    println!("number of tests = {}", NUMBER_OF_TESTS);
    println!("sieve error margin = {}", SIEVE_ERROR_MARGIN);

    println!();

    println!("right guess rate = {:.1}", number_of_right_guess as f32 / NUMBER_OF_TESTS as f32 * 100.0);
    println!("wrong guess rate = {:.1}", number_of_wrong_guess as f32 / NUMBER_OF_TESTS as f32 * 100.0);
    println!("rate of errors (can't find enough B-smooth) = {:.1}", number_of_error_cant_find_enough_b_smooth as f32 / NUMBER_OF_TESTS as f32 * 100.0);
    println!("average number of trails to find zero vector = {:.1}", number_of_trails_to_find_zero_vector as f32 / NUMBER_OF_TESTS as f32 * 100.0);
    
    println!();

    println!("rate of b-smooth numbers found out of all tested numbers = {:.3}", number_of_actual_b_smooth_numbers as f32 / number_of_tested_numbers as f32 * 100.0);
    println!("rate of b-smooth numbers found out of sieved numbers = {:.3}", number_of_actual_b_smooth_numbers as f32 / number_of_sieved_numbers as f32 * 100.0);
    
    println!();

    println!("average computation time = {:#?}", average_time_of_function)
}

