
def numeric_user_input(lower_range: int, upper_range: int, first_prompt: str,
                       following_prompt="") -> int:
    """
    Very basic input checking interface where program asks for a number from a
    range of contiguous ones.

    :param lower_range: int
    :param upper_range: int
    :param first_prompt: str
    :param following_prompt: only used if invalid user input was entered
    :return:
    """
    user_input = int(input(first_prompt))
    while user_input < lower_range or user_input > upper_range:
        if following_prompt == "":
            user_input = int(input("Select between " + str(lower_range) +
                                   " and " + str(upper_range) + ", inclusive."))
        else:
            user_input = int(input(following_prompt))

    return user_input


