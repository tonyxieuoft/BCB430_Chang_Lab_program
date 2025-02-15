

def expect_threshold_user_input() -> float:
    """
    Same as numeric_user_input, but for expect threshold values
    :return: the user-specified expect threshold
    """
    user_input = float(input())
    while user_input <= 0:
        user_input = float(input("Incorrect. Please enter a value greater"
                                 "than 0 and less than or equal to 1"))

    return user_input
