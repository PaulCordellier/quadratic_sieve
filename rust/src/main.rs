use std::{fs::File, io::Write, time::{Duration, Instant}};

mod quadratic_sieve;
mod sieve_of_eratosthenes;
mod section_2_3;
mod get_zero_vector_combination;

fn main() {

    let content = "exponent of n;\
            computation time;\
            rate of successful results;\
            average number of trails to find zero vector;\
            rate of sieved numbers that are B-smooth";

    let mut file = File::create("data.csv").expect("Can not create or override data.csv");
    
    file.write_all(content.as_bytes()).expect("Can not write on the file");

    for i in 8..29 {

        println!("testing for exponent {}", i);

        let benckmark_results = benckmark_quadratic_sieve(i);

        let csv_line: String = format!(
            "\n{};{:.3};{:.4};{};{:.4}",
            i,
            benckmark_results.average_computation_time.as_secs_f64() * 100.0,
            benckmark_results.rate_of_successful_results,
            benckmark_results.average_number_of_trails_to_find_zero_vector,
            if benckmark_results.rate_of_sieved_numbers_that_are_b_smooth.is_nan() { 
                0.0 
            } else {
                benckmark_results.rate_of_sieved_numbers_that_are_b_smooth 
            }
        );

        file.write(csv_line.as_bytes()).expect("Can not write on the file");
    }
}

struct BenckmarkResults {
    average_computation_time: Duration,
    rate_of_successful_results: f32,
    average_number_of_trails_to_find_zero_vector: u64,
    rate_of_sieved_numbers_that_are_b_smooth: f32
}

fn benckmark_quadratic_sieve(exponent_of_n: u32) -> BenckmarkResults {

    let mut benckmark_results = BenckmarkResults {
        average_computation_time: Duration::new(0, 0),
        rate_of_successful_results: 0.0,
        average_number_of_trails_to_find_zero_vector: 0,
        rate_of_sieved_numbers_that_are_b_smooth: 0.0
    };

    let test_start_number: u128 = 10_u128.pow(exponent_of_n) + 1;
    const NUMBER_OF_TESTS: u128 = 50;
    const SIEVE_MARGIN_OF_ERROR: u32 = 20;

    let mut number_of_actual_b_smooth_numbers = 0;
    let mut number_of_sieved_numbers = 0;

    for n in (test_start_number .. test_start_number + NUMBER_OF_TESTS * 2).step_by(2) {

        let time_before_quadratic_sieve = Instant::now();

        let algo_results = quadratic_sieve::quadratic_sieve(
            n, 
            SIEVE_MARGIN_OF_ERROR, 
            3
        );

        let elapsed_time = time_before_quadratic_sieve.elapsed();

        let algo_results = match algo_results {
            Ok(algo_results) => algo_results,
            Err(_error_message) => {
                continue;
            }
        };

        if algo_results.non_trivial_factor == 1 {
            continue;
        }
        if n % algo_results.non_trivial_factor != 0 {
            panic!("There is a problem in the quadratic sieve algorithm: {} isn't a factor of {}", algo_results.non_trivial_factor, n);
        }
        
        benckmark_results.average_computation_time += elapsed_time;

        benckmark_results.rate_of_successful_results += 1.0;

        number_of_actual_b_smooth_numbers += algo_results.nb_of_b_smooth_found;
        number_of_sieved_numbers += algo_results.nb_of_sieved_numbers;
        benckmark_results.average_number_of_trails_to_find_zero_vector += algo_results.nb_of_trails_to_find_zero_vector;
    }

    benckmark_results.average_computation_time /= NUMBER_OF_TESTS as u32;
    benckmark_results.average_number_of_trails_to_find_zero_vector /= NUMBER_OF_TESTS as u64;
    benckmark_results.rate_of_sieved_numbers_that_are_b_smooth = number_of_actual_b_smooth_numbers as f32 / number_of_sieved_numbers as f32;

    benckmark_results
}

