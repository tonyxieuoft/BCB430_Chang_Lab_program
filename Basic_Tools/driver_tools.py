import time
from typing import List

from selenium.common import WebDriverException
from selenium.webdriver.chrome import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.remote.webelement import WebElement
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

WAIT_CONSTANT = 0.01


class DriverTimeoutException(Exception):
    pass


def get_element(specific_driver, by_what, description, timer_limit) -> \
        WebElement:
    """
    Wait until an element is present, then get the element as an WebElement.

    :param specific_driver: driver, which allows for narrowing of the scope of
    elements searched for
    :param by_what: By.ID, By.XPATH, etc. the type of the descriptor used to
    search for the element
    :param description: the description of the element, ex. its xpath
    :param timer_limit: amount of time to attempt for
    :return: the element as an WebElement
    """

    counter = 0
    while counter < timer_limit:
        try:
            element = specific_driver.find_element(by_what, description)
            return element
        except WebDriverException:
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    # throws an exception if nothing is got within timer limit
    raise DriverTimeoutException(
        "Timeout error for getting element, description: " +
        description)


def get_elements(driver: webdriver, by_what: str, description: str,
                 timer_limit: int) -> List:
    WebDriverWait(driver, timer_limit).until(
        EC.presence_of_element_located((by_what, description)))
    elements = driver.find_elements(by_what, description)

    return elements


def try_click(specific_driver, by_what, description) -> bool:
    """
    Attempts a click, returning True if successful and False otherwise.

    :param specific_driver: driver
    :param by_what: type of description, ex. By.ID, By.XPATH
    :param description: description of element to click, in string form
    :return: True iff the click was successful
    """

    try:
        specific_driver.find_element(by_what, description).click()
        return True
    except WebDriverException:
        return False


def clicking_wrapper(specific_driver, by_what, description,
                     timer_limit) -> None:
    """
    Keeps attempting to click an element until timer_limit is reached or
    success occurs.

    :param specific_driver: a more specific driver that limits the scope of
    elements considered
    :param by_what: the type of the description (ex. By.ID, By.XPATH)
    :param description: descriptor of the element to click, as a string
    :param timer_limit: stop attempting once this much time has passed.
    :return:
    """

    counter = 0
    while counter < timer_limit:
        try:
            # finds the element and attempts to click
            specific_driver.find_element(by_what, description).click()
            break
        except WebDriverException:
            # pauses for a period of time
            time.sleep(WAIT_CONSTANT)
            counter += WAIT_CONSTANT

    if counter >= timer_limit:
        # throw an exception if time limit is reached
        raise DriverTimeoutException(
            "Timeout error for clicking wrapper, description:" +
            description)
