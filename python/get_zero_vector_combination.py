import copy;

vector_addition_info = tuple[list[int], list[int]]  # Creation type used in the function zero_vector_is_found

# Function for step 3:
def get_zero_vector_combination(vector_matrix_mod_2: list[list[int]], K: int) -> list[int] | None:

    # The type vector_addition_info is a tuple composed of the indexes of the vectors that were added,
    # and the result of the additions.

    vector_infos: list[vector_addition_info] = [([i], copy.deepcopy(vector_matrix_mod_2[i])) for i in range(len(vector_matrix_mod_2))]

    for i in range(K):

        vector_info_to_add: vector_addition_info = ()

        # looking for addition results of vectors where addition_result[i] == 1:

        for j in range(len(vector_infos)):
            if vector_infos[j][1][i] == 1:
                vector_info_to_add = vector_infos.pop(j)
                break

        if vector_info_to_add == ():
            continue

        for j in range(len(vector_infos)):

            if vector_infos[j][1][i] == 0:
                continue

            append_lists_without_duplicates(vector_infos[j][0], vector_info_to_add[0])

            vector_addition = [0 if vector_infos[j][1][k] == vector_info_to_add[1][k] else 1 for k in range(len(vector_info_to_add[1]))]

            vector_infos[j] = (vector_infos[j][0], vector_addition)

    indexes_for_zero_vector = []

    for (indexes_of_added_vectors, added_vector) in vector_infos:

        if all(value == 0 for value in added_vector):
            append_lists_without_duplicates(indexes_for_zero_vector, indexes_of_added_vectors)

    if len(indexes_for_zero_vector) >= K + 1:
        return indexes_for_zero_vector

    return None


def append_lists_without_duplicates(list_with_added_elements: list, list_to_add: list):
    """
    If the list has the same element two times, this element won't appear in the
    list "list_with_added_elements".
    """
    
    for i in list_to_add:
        if i in list_with_added_elements:
            list_with_added_elements.remove(i)
        else:
            list_with_added_elements.append(i)

    # this other version of this function is somehow slower: its aim was to cancel
    # the append method if it removes more element than it adds

    # elements_to_add = []
    # elements_to_remove = []

    # for i in list_to_add:
    #     if i in list_with_added_elements:
    #         elements_to_remove.append(i)
    #     else:
    #         elements_to_add.append(i)

    # if len(elements_to_add) >= len(elements_to_remove):
    #     list_with_added_elements.extend(elements_to_add)
    #     for i in elements_to_remove:
    #         list_with_added157289988114337_elements.remove(i)