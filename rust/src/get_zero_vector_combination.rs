use std::vec;

/// The goal is to find a combination for a zero vector form a matrix of vectors mod 2. Returns a
/// vector of indexes, or nothing. When the vector of indexes is returned, the indexes represent
/// the vectors that, added, will result in a zero vector.
#[allow(non_snake_case)]
pub fn get_zero_vector_combination(vector_matrix_mod_2: &Vec<Vec<u32>>, K: usize) -> Option<Vec<usize>> {

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

    // Here we perform the gaussian elimination:

    for i in 0..K + 1 {
        
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

    // At this point, every addition result is zero vector.
    // We need here to find the longest combination of indexes that added result in the zero vector

    let mut matrix_of_indexes_for_zero_vector: Vec<Vec<usize>> = vec![vec![]];

    for vector_info in vector_infos {
        
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

    assert_eq!(
        vec_with_added_elements.len(),
        vec_to_add.len(),
        "vector_addition_mod_2: The length of the two vectors should be the same to make the vector addition possible."
    );
    
    for i in 0..vec_to_add.len() {
        vec_with_added_elements[i] = if vec_with_added_elements[i] == vec_to_add[i] { 0 } else { 1 };
    }
}