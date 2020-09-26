import numpy as np
import pandas as pd


# header_flag = 1 & index_flag = 1 -> rows and columns both labeled
# header_flag = None & index_flag = None -> rows and columns both unlabeled
def load_file(path, header_flag, index_flag):
    """
    :param path: Path to file
    :param header_flag:
    :param index_flag:
    :return: Loaded matrix or list as pandas dataframe
    """
    file = pd.read_csv(path, header=header_flag, index_col=index_flag)
    return file

def sum_col(matrix):
    """
    :param matrix: A matrix
    :return: Array of column sums of input matrix
    """
    return np.array([sum(tf) for tf in zip(*matrix)])
