import unittest

import numpy as np

from src.util.train_test_split import GroupShuffleSplitND


class TestGroupShuffleSplitND(unittest.TestCase):
    def setUp(self) -> None:
        self.X = np.arange(10000, dtype="int")
        self.groups_2d = np.random.randint(10, size=(10000, 2), dtype="int")
        self.groups_3d = np.random.randint(10, size=(10000, 3), dtype="int")
        self.splitter = GroupShuffleSplitND(
            n_splits=5, test_size=0.001, random_state=np.random.RandomState(42)
        )

    def test_no_overlap_between_train_and_test_for_2d(self):
        """Check that no samples are in both train and test set for 2D groups."""
        for i, (train, test) in enumerate(
            self.splitter.split(self.X, groups=self.groups_2d)
        ):
            with self.subTest(fold=i):
                self.assertFalse(set(train).intersection(set(test)))

    def test_test_groups_not_in_train_for_2d(self):
        """Check that no groups are in both train and test set for 2D groups."""
        for i, (train, test) in enumerate(
            self.splitter.split(self.X, groups=self.groups_2d)
        ):
            test_groups = self.groups_2d[test, :]
            train_groups = self.groups_2d[train, :]
            for test_group in test_groups:
                with self.subTest(test_group=test_group, fold=i):
                    self.assertFalse((train_groups == test_group).all(axis=1).any())

    def test_no_overlap_between_train_and_test_for_3d(self):
        """Check that no samples are in both train and test set for 3D groups."""
        for i, (train, test) in enumerate(
            self.splitter.split(self.X, groups=self.groups_3d)
        ):
            with self.subTest(fold=i):
                self.assertFalse(set(train).intersection(set(test)))

    def test_test_groups_not_in_train_for_3d(self):
        """Check that no groups are in both train and test set for 3D groups."""
        for i, (train, test) in enumerate(
            self.splitter.split(self.X, groups=self.groups_3d)
        ):
            test_groups = self.groups_3d[test, :]
            train_groups = self.groups_3d[train, :]
            for test_group in test_groups:
                with self.subTest(test_group=test_group, fold=i):
                    self.assertFalse((train_groups == test_group).all(axis=1).any())


if __name__ == "__main__":
    unittest.main()
