from functools import reduce

import numpy as np
from sklearn.model_selection import ShuffleSplit
from sklearn.utils.validation import check_array


class GroupShuffleSplitND(ShuffleSplit):
    """
    Extension of scikit-learn's GroupShuffleSplit to apply to n-dimensional groups.

    Splits data which is grouped along n dimensions (e.g. data from n-component chemical reactions) so that the
    train and test set groups do not overlap on any dimension. Note that this typically means that some data will not
    be used, neither in the train, nor test set. F
    The train_size or test_size arguments refer to the ratio of groups and not the ratio of samples.
    The ratio is applied along each dimension separately, so if the test_size is 0.1, and n = 3, then the expected
    size of the test set is 0.1^3 = 0.1% of the total data, if all groups are of equal size.
    """

    def __init__(
        self, n_splits=5, *, test_size=None, train_size=None, random_state=None
    ):
        super().__init__(
            n_splits=n_splits,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state,
        )
        self._default_test_size = 0.1

    def _iter_indices(self, X, y, groups):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        # validate input
        groups = check_array(groups, ensure_2d=True, dtype=None)

        # determine the unique group labels present along each dimension
        n_groups = groups.shape[1]
        group_labels, group_indices = zip(
            *[(np.unique(groups[:, i], return_inverse=True)) for i in range(n_groups)]
        )

        # split the groups into train and test set for each dimension
        # note: we need to use the 2-argument form of super() here, bc list comprehension has its own nested scope
        folds = list(
            zip(
                *[
                    super(GroupShuffleSplitND, self)._iter_indices(
                        X=list(range(len(group_labels[i])))
                    )
                    for i in range(n_groups)
                ]
            )
        )

        for fold in folds:
            group_ind_train, group_ind_test = zip(*fold)
            # determine datapoints that are in the test and train set for each dimension
            test_ind = [
                X[np.isin(group_ind, ind_test)]
                for group_ind, ind_test in zip(group_indices, group_ind_test)
            ]
            train_ind = [
                X[np.isin(group_ind, ind_train)]
                for group_ind, ind_train in zip(group_indices, group_ind_train)
            ]

            # aggregate over all dimensions:
            # only indices that appear in a set for all dimensions are included in the final set
            test = reduce(np.intersect1d, test_ind)
            train = reduce(np.intersect1d, train_ind)

            # check that we are not returning empty sets
            if len(train) == 0:
                raise RuntimeError(
                    "No samples found in train groups. Consider lowering test_size."
                )
            if len(test) == 0:
                raise RuntimeError(
                    "No samples found in test groups. This can occur because test_size is to small, or, due to chance if some group combination contain zero samples. Try increasing test_size or use a different random seed."
                )
            yield train, test

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : array-like of shape (n_samples,), default=None
            The target variable for supervised learning problems.

        groups : array-like of shape (n_samples, n_groups)
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.

        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting `random_state`.
        Note that the random_state is used n_groups times.
        If n_groups > 1 and random_state is an integer, then the split for all groups use the same seed.
        To avoid this, pass a numpy RandomState object.
        """
        return super().split(X, y, groups)
