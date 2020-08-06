import os
import time
from multiprocessing import Process


def get_temp_directory():
    """
    create temp directory
    """
    current = os.getcwd()
    temp = os.path.join(current, "temp")
    os.makedirs(temp, exist_ok=True)
    return temp


def remove_temp_data(temp):
    """
    Remove the data in the temp directory (param: temp)
    """
    walker = os.walk(temp, topdown=False)
    file_names = []
    dir_names = []
    for root, dirs, files in walker:
        for file_ in files:
            name = os.path.join(root, file_)
            file_names.append(name)
        for dir_ in dirs:
            name = os.path.join(root, dir_)
            dir_names.append(name)
    # print(file_names)
    for file_ in file_names:
        os.remove(file_)
    for dir_ in dir_names:
        os.rmdir(dir_)
    os.rmdir(temp)
    

def timeout_function(func, kw, timeout):
    """
    Submit a function and its kwargs to a dask client. Wait for a
    number of seconds before cancelling the future and closing the
    client.

    Parameters
    ----------
    func: function
    kw: dict
        key word arguments required to run said function
    timeout: int
        number of seconds for which to wait before cancelling the future
    """
    f = func.__name__
    p = Process(target=func, kwargs=kw, name=f)
    p.start()
    time.sleep(timeout)
    p.terminate()
    if not p.is_alive():
        m = f'The function {f} was cancelled after {timeout} seconds'
        print(m)
