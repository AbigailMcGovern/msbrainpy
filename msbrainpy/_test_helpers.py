import os
import time
from dask.distributed import Client

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
    files = os.listdir(temp)
    for file_ in files:
        path = os.path.join(temp, file_)
        os.remove(path)
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
    c = Client()
    future = c.submit(func, **kw)
    time.sleep(timeout)
    c.cancel(future)
    f = func.__name__
    if future.cancelled():
        m = f'The function {f} was cancelled after {timeout} seconds'
        print(m)
    c.close()