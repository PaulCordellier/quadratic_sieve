use std::vec;

#[allow(non_snake_case)]
pub fn get_zero_vector_combination(vector_matrix_mod_2: &Vec<Vec<u32>>, K: usize) -> Option<Vec<usize>> {

    #[derive(Debug)]
    struct VectorInfo {
        indices_of_added_vectors: Vec<usize>,
        addition_result: Vec<u32>
    }

    let enumerated_vectors = vector_matrix_mod_2.iter().enumerate();
    // enumerated_vectors is an iterator of the vectors mod 2 with their indexes
        
    let mut vector_infos: Vec<VectorInfo> = enumerated_vectors.map(|(index, vector_mod_2)| {
        VectorInfo {
            indices_of_added_vectors: vec![index],
            addition_result: vector_mod_2.clone()
        }
    })
    .collect();

    for i in 0..K {
        
        let mut vector_info_to_add: Option<VectorInfo> = None;

        // looking for addition results of vectors where addition_result[i] == 1:

        for j in 0..vector_infos.len() {
            if vector_infos[j].addition_result[i] == 1 {
                vector_info_to_add = Some(vector_infos.swap_remove(j));
                break;
            }
        }

        let vector_info_to_add: VectorInfo = match vector_info_to_add {
            Some(x) => x,
            None => continue
        };

        for vector_info in vector_infos.iter_mut() {

            if vector_info.addition_result[i] == 0 {
                continue;
            }

            append_index_vecs_without_duplicates(
                &mut vector_info.indices_of_added_vectors,
                &vector_info_to_add.indices_of_added_vectors
            );

            vector_addition_mod_2(
                &mut vector_info.addition_result,
                &vector_info_to_add.addition_result
            );
        }
    }

    let mut matrix_of_indexes_for_zero_vector: Vec<Vec<usize>> = vec![vec![]];

    for vector_info in vector_infos {
        
        if vector_info.addition_result.iter().any(|&x| x != 0) {
            continue;
        }

        for i in 0..matrix_of_indexes_for_zero_vector.len() {
            
            let mut vector_of_indexes = matrix_of_indexes_for_zero_vector[i].clone();
        
            append_index_vecs_without_duplicates(
                &mut vector_of_indexes,
                &vector_info.indices_of_added_vectors
            );

            matrix_of_indexes_for_zero_vector.push(vector_of_indexes);
        }

        // Keep only the file longest vectors in the list:

        const MATRIX_MAX_SIZE: usize = 30;

        if matrix_of_indexes_for_zero_vector.len() <= MATRIX_MAX_SIZE {
            continue;
        }

        matrix_of_indexes_for_zero_vector.sort_by(|vec1, vec2| vec2.len().cmp(&vec1.len()));
        matrix_of_indexes_for_zero_vector.truncate(MATRIX_MAX_SIZE);
    }

    let longest_vector_of_indices = matrix_of_indexes_for_zero_vector.into_iter()
                                                .next()
                                                .expect("matrix_of_indexes_for_zero_vector shouldn't be empty here.");

    if longest_vector_of_indices.len() >= K + 1 {
        return Some(longest_vector_of_indices);
    }

    None
}

fn append_index_vecs_without_duplicates(vec_with_added_elements: &mut Vec<usize>, vec_to_add: &Vec<usize>) {
    for &i in vec_to_add {

        let position_of_i: Option<usize> = vec_with_added_elements.iter().position(|&x| x == i);

        match position_of_i {
            Some(position_of_i) => {
                vec_with_added_elements.swap_remove(position_of_i);
            },
            None => vec_with_added_elements.push(i)
        }
    }
}

fn vector_addition_mod_2(vec_with_added_elements: &mut Vec<u32>, vec_to_add: &Vec<u32>) {

    let vector_len = vec_with_added_elements.len();
    
    if vector_len != vec_to_add.len() {
        panic!("The length of the two vectors should be the same to make the vector addition possible.")
    }

    for i in 0..vector_len {
        vec_with_added_elements[i] = if vec_with_added_elements[i] == vec_to_add[i] { 0 } else { 1 };
    }
}