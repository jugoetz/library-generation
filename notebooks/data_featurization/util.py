import numpy as np
import pandas as pd

from src.definitions import DATA_DIR


def write_indices_and_stats(
    indices,
    sizes,
    pos_class,
    split_dimension,
    train_size,
    total_size,
    data_name,
    n_initiators=None,
    n_monomers=None,
    n_terminators=None,
    save_indices=True,
):
    """
    Write function that can be reused for the other splits
    Args:
        indices: list of 3-tuples
        sizes: list of 3-tuples, length equal to indices
        pos_class: list of 3-tuples, length equal to indices, containing three Sequences with length == n_labels
        split_dimension: str or int, e.g. "0", "1", "2", "3"
        train_size: ratio of train samples from all samples. Give as percentage. If _train_size suffix should be omitted, pass an empty string
        total_size (int): number of samples in the data set
        data_name (str): name of the dataset. Used as base for the path where indices and stats will be saved.
        n_initiators: list of 3-tuples of int, number of distinct initiators
        n_monomers: same as n_initiators, for monomers
        n_terminators: same as n_initiators, for terminators
        save_indices: bool, whether to save the indices. Useful if we need to regenerate statistics. Default: True

    Returns:
        None
    """
    n_folds = len(indices)

    save_dir = (
        DATA_DIR
        / "curated_data"
        / "splits"
        / f"{data_name}_{split_dimension}D_split_{train_size}".rstrip("_")
    )
    save_dir.mkdir(parents=True, exist_ok=True)
    if save_indices:
        for i, (idx_train, idx_val, idx_test) in enumerate(indices):
            with open(save_dir / f"fold{i}_train.csv", "w") as f:
                f.write("index\n")
                f.write("\n".join([str(i) for i in idx_train]))

            with open(save_dir / f"fold{i}_val.csv", "w") as f:
                f.write("index\n")
                f.write("\n".join([str(i) for i in idx_val]))

            with open(save_dir / f"fold{i}_test.csv", "w") as f:
                f.write("index\n")
                f.write("\n".join([str(i) for i in idx_test]))

    for i, (size, count_pos) in enumerate(zip(sizes, pos_class)):
        with open(save_dir / f"fold{i}_statistics.txt", "w") as f:
            f.write(f"Train samples: {size[0]} ({size[0]/total_size:.1%})\n")
            f.write(f"Val samples: {size[1]} ({size[1]/total_size:.1%})\n")
            f.write(f"Test samples: {size[2]} ({size[2]/total_size:.1%})\n")
            if split_dimension > 1:
                f.write(
                    f"Not used: {total_size - np.sum(size):.0f} ({(total_size - np.sum(size)) / total_size:.1%})\n"
                )
            f.write(
                f"Train samples binary_A has label 1: {count_pos[0][0]} ({count_pos[0][0]/size[0]:.1%})\n"
            )
            f.write(
                f"Train samples binary_B has label 1: {count_pos[0][1]} ({count_pos[0][1]/size[0]:.1%})\n"
            )
            f.write(
                f"Train samples binary_C has label 1: {count_pos[0][2]} ({count_pos[0][2]/size[0]:.1%})\n"
            )

            f.write(
                f"Val samples binary_A has label 1: {count_pos[1][0]} ({count_pos[1][0]/size[1]:.1%})\n"
            )
            f.write(
                f"Val samples binary_B has label 1: {count_pos[1][1]} ({count_pos[1][1]/size[1]:.1%})\n"
            )
            f.write(
                f"Val samples binary_C has label 1: {count_pos[1][2]} ({count_pos[1][2]/size[1]:.1%})\n"
            )

            f.write(
                f"Test samples binary_A has label 1: {count_pos[2][0]} ({count_pos[2][0]/size[2]:.1%})\n"
            )
            f.write(
                f"Test samples binary_B has label 1: {count_pos[2][1]} ({count_pos[2][1]/size[2]:.1%})\n"
            )
            f.write(
                f"Test samples binary_C has label 1: {count_pos[2][2]} ({count_pos[2][2]/size[2]:.1%})\n"
            )

            f.write(
                f"Chance level average precision macro on val set: {np.sum(count_pos[1]) / (3 * size[1]):.3f}\n"
            )
            f.write(
                f"Chance level average precision macro on test set: {np.sum(count_pos[2]) / (3 * size[2]):.3f}\n"
            )

    if n_initiators:
        for i, n in enumerate(n_initiators):
            with open(save_dir / f"fold{i}_statistics.txt", "a") as f:
                f.write(f"Train initiators: {n[0]}\n")
                f.write(f"Val initiators: {n[1]}\n")
                f.write(f"Test initiators: {n[2]}\n")

    if n_monomers:
        for i, n in enumerate(n_monomers):
            with open(save_dir / f"fold{i}_statistics.txt", "a") as f:
                f.write(f"Train monomers: {n[0]}\n")
                f.write(f"Val monomers: {n[1]}\n")
                f.write(f"Test monomers: {n[2]}\n")

    if n_terminators:
        for i, n in enumerate(n_terminators):
            with open(save_dir / f"fold{i}_statistics.txt", "a") as f:
                f.write(f"Train terminators: {n[0]}\n")
                f.write(f"Val terminators: {n[1]}\n")
                f.write(f"Test terminators: {n[2]}\n")

    # summary statistics
    sum_pos_class = np.sum(pos_class, axis=0)
    sum_sizes = np.sum(sizes, axis=0)
    with open(save_dir / "summary_statistics.txt", "w") as f:
        f.write(
            f"Mean Train samples: {sum_sizes[0] / n_folds:.0f} ({sum_sizes[0] / n_folds / total_size:.1%})\n"
        )
        f.write(
            f"Mean Val samples: {sum_sizes[1] / n_folds:.0f} ({sum_sizes[1] / n_folds / total_size:.1%})\n"
        )
        f.write(
            f"Mean Test samples: {sum_sizes[2] / n_folds:.0f} ({sum_sizes[2] / n_folds / total_size:.1%})\n"
        )
        if split_dimension > 1:
            f.write(
                f"Not used: {total_size - np.sum(sum_sizes) / n_folds:.0f} ({(total_size - np.sum(sum_sizes) / n_folds) / total_size:.1%})\n"
            )
        f.write(
            f"Mean Train samples binary_A has label 1: {sum_pos_class[0][0] / n_folds:.0f} ({sum_pos_class[0][0]/sum_sizes[0]:.1%})\n"
        )
        f.write(
            f"Mean Train samples binary_B has label 1: {sum_pos_class[0][1] / n_folds:.0f} ({sum_pos_class[0][1]/sum_sizes[0]:.1%})\n"
        )
        f.write(
            f"Mean Train samples binary_C has label 1: {sum_pos_class[0][2] / n_folds:.0f} ({sum_pos_class[0][2]/sum_sizes[0]:.1%})\n"
        )

        f.write(
            f"Mean Val samples binary_A has label 1: {sum_pos_class[1][0] / n_folds:.0f} ({sum_pos_class[1][0]/sum_sizes[1]:.1%})\n"
        )
        f.write(
            f"Mean Val samples binary_B has label 1: {sum_pos_class[1][1] / n_folds:.0f} ({sum_pos_class[1][1]/sum_sizes[1]:.1%})\n"
        )
        f.write(
            f"Mean Val samples binary_C has label 1: {sum_pos_class[1][2] / n_folds:.0f} ({sum_pos_class[1][2]/sum_sizes[1]:.1%})\n"
        )

        f.write(
            f"Mean Test samples binary_A has label 1: {sum_pos_class[2][0] / n_folds:.0f} ({sum_pos_class[2][0]/sum_sizes[2]:.1%})\n"
        )
        f.write(
            f"Mean Test samples binary_B has label 1: {sum_pos_class[2][1] / n_folds:.0f} ({sum_pos_class[2][1]/sum_sizes[2]:.1%})\n"
        )
        f.write(
            f"Mean Test samples binary_C has label 1: {sum_pos_class[2][2] / n_folds:.0f} ({sum_pos_class[2][2]/sum_sizes[2]:.1%})\n"
        )

        f.write(
            f"Mean chance level average precision macro on val set: {np.sum(sum_pos_class[1]) / (3 * sum_sizes[1]):.3f}\n"
        )
        f.write(
            f"Mean chance level average precision macro on test set: {np.sum(sum_pos_class[2]) / (3 * sum_sizes[2]):.3f}\n"
        )

        if n_initiators:
            mean = np.mean(n_initiators, axis=0)
            f.write(f"Mean Train initiators: {mean[0]}\n")
            f.write(f"Mean Val initiators: {mean[1]}\n")
            f.write(f"Mean Test initiators: {mean[2]}\n")
        if n_monomers:
            mean = np.mean(n_monomers, axis=0)
            f.write(f"Mean Train monomers: {mean[0]}\n")
            f.write(f"Mean Val monomers: {mean[1]}\n")
            f.write(f"Mean Test monomers: {mean[2]}\n")
        if n_terminators:
            mean = np.mean(n_terminators, axis=0)
            f.write(f"Mean Train terminators: {mean[0]}\n")
            f.write(f"Mean Val terminators: {mean[1]}\n")
            f.write(f"Mean Test terminators: {mean[2]}\n")
