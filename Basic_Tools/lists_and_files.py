import os
import time
from typing import Iterable, List, Dict, Tuple
from Bio import Entrez
import re


def file_to_list(filepath: str) -> List:
    """
    Converts the contents of a file into a list by line

    :param filepath: the path of the file to convert
    :return: a list with each index storing a line of the file
    """

    f = open(filepath)
    raw = f.readlines()
    output = []
    for i in range(0, len(raw)):

        # strips to remove any potential whitespace
        if raw[i].strip() != "":
            output.append(raw[i].strip())

    return output


def list_to_string(lst: Iterable[object], delim) -> str:
    """
    Converts a list to a string, with each index separated by a comma

    :param lst: the list to convert
    :return: a comma-separated string
    """
    first = True
    string = ""
    for item in lst:
        if first:
            string = str(item)
            first = False
        else:
            string = string + delim + str(item)

    return string


def unique_filepath(file_path: str) -> str:
    """
    Return a unique filename by adding an increment version to the end of the
    provided filename. For example, if the file 'filename' exists, this function
    will check to see if 'filename (1)' exists as a name and return it if not.
    Otherwise, will keep incrementing the version to 'filename (2)',
    'filename (3)', etc.

    :param file_path: name of the file to find a unique name for
    :return: the filename with a unique increment version attached to the end
    """
    if not os.path.exists(file_path):
        return file_path

    else:
        path_sections = os.path.splitext(file_path)
        iteration = 1
        altered_name = path_sections[0] + str(iteration) + path_sections[1]
        while os.path.exists(altered_name):
            iteration += 1
            altered_name = path_sections[0] + str(iteration) + path_sections[1]

        return altered_name


def make_unique_directory(parent, name):

    dir_path = os.path.join(parent, name)
    dir_path = unique_filepath(dir_path)
    os.mkdir(dir_path)

    return dir_path
