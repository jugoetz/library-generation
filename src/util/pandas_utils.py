import numpy as np
import pandas as pd


def alphanumeric_index_sort(index):
    """custom sorting function to sort an alphanumeric index with one letter and arbitrary digits only by the number part"""
    return index.str[1:].astype(int)


def sort_by_sparsity(df: pd.DataFrame, axis: int = 0) -> pd.DataFrame:
    """
    Sorts a dataframe with NaN values by sparsity of its columns or rows.
    The resulting dataframe will have NaN entries clustered on one side.

    Note: we do not sort by the number of NaNs, but by their position.

    By applying this function twice with different values for the `axis` parameter, a pivoted dataframe or table can be
    sorted by sparsity.

    Args:
        df (pd.DataFrame): The dataframe to be sorted.
        axis (int): The axis along which to sort. Must be 0 or 1.

    Returns:
        pd.DataFrame: The sorted dataframe.

    Raises:
        ValueError: If the axis parameter is not 0 or 1.
    """
    if axis == 0:
        df_tmp = df.copy()
    elif axis == 1:
        df_tmp = df.T.copy()
    else:
        raise ValueError
    # here we will store rows for the new df
    rows = []
    while len(df_tmp) > 0:
        # we take the first row and identify all columns that are not nan
        columns = df_tmp.columns[~df_tmp.iloc[0, :].isna().values]
        # we identify all rows that have the same pattern
        rows.append(df_tmp[columns].dropna(axis=0, how="any"))
        # we delete all just collected rows from df. Errors are ignored because from the second iteration on, this will
        # try to delete rows that are already deleted.
        df_tmp = df_tmp.drop(
            np.concatenate([r.index.values for r in rows]), errors="ignore"
        )
    if axis == 0:
        return pd.concat(rows)
    elif axis == 1:
        return pd.concat(rows).T
