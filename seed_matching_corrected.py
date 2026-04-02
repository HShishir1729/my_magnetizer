import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import linear_sum_assignment
import os

base = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new"
galaxy_file = "/home/shishir/test/Magnetizer/params_template/galaxy_list10.txt"


def read_variable_rows(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.strip():
                data.append(np.array(line.split(), dtype=float))
    return data


galaxies = np.loadtxt(galaxy_file, dtype=int)

print("Total galaxies:", len(galaxies))

for gal_index in galaxies:

    print("\n==============================")
    print("Processing galaxy:", gal_index)
    print("==============================")

    seed_bp_file = f"{base}/seed/seed_output_{gal_index}_fid_floor.dat"
    seed_r_file  = f"{base}/r_seed/radius_seed_{gal_index}_fid_floor.dat"

    raw_bp_file = f"{base}/bp/Bp_profiles_{gal_index}_fid_floor.dat"
    raw_r_file  = f"{base}/r/r_profiles_{gal_index}_fid_floor.dat"

    if not (
        os.path.exists(seed_bp_file)
        and os.path.exists(seed_r_file)
        and os.path.exists(raw_bp_file)
        and os.path.exists(raw_r_file)
    ):
        print("Missing files. Skipping galaxy.")
        continue

    # -------------------------
    # READ SEED FILES
    # -------------------------
    Bp_profiles_seed = read_variable_rows(seed_bp_file)
    r_profiles_seed  = read_variable_rows(seed_r_file)

    nrows_seed = len(r_profiles_seed)

    Bp_sampled_seed = np.zeros((nrows_seed, 5))
    r_sampled_seed  = np.zeros((nrows_seed, 5))

    for i in range(nrows_seed):

        r = r_profiles_seed[i]
        Bp = Bp_profiles_seed[i]

        zero_idx = np.where(np.isclose(Bp, 0))[0]

        if len(zero_idx) < 2:
            continue

        start = zero_idx[0]
        end   = zero_idx[-1]

        r_valid = r[start:end+1]
        Bp_valid = Bp[start:end+1]

        sort_idx = np.argsort(r_valid)

        r_valid = r_valid[sort_idx]
        Bp_valid = Bp_valid[sort_idx]

        r_disk = np.max(r_valid) / 2.7
        r_eval = r_disk * np.array([1/4, 2/4, 3/4, 4/4, 5/4])

        interp_func = interp1d(
            r_valid,
            Bp_valid,
            kind="linear",
            fill_value="extrapolate"
        )

        Bp_sampled_seed[i] = interp_func(r_eval)
        r_sampled_seed[i]  = r_eval

    # -------------------------
    # READ RAW FILES
    # -------------------------
    Bp_profiles_raw = read_variable_rows(raw_bp_file)
    r_profiles_raw  = read_variable_rows(raw_r_file)

    nrows_raw = len(r_profiles_raw)

    Bp_sampled_raw = np.zeros((nrows_raw, 5))
    r_sampled_raw  = np.zeros((nrows_raw, 5))

    for i in range(nrows_raw):

        r = r_profiles_raw[i]
        Bp = Bp_profiles_raw[i]

        zero_idx = np.where(np.isclose(Bp, 0))[0]

        if len(zero_idx) < 2:
            continue

        start = zero_idx[0]
        end   = zero_idx[-1]

        r_valid = r[start:end+1]
        Bp_valid = Bp[start:end+1]

        sort_idx = np.argsort(r_valid)

        r_valid = r_valid[sort_idx]
        Bp_valid = Bp_valid[sort_idx]

        r_disk = np.max(r_valid) / 2.7
        r_eval = r_disk * np.array([1/4, 2/4, 3/4, 4/4, 5/4])

        interp_func = interp1d(
            r_valid,
            Bp_valid,
            kind="linear",
            fill_value="extrapolate"
        )

        Bp_sampled_raw[i] = interp_func(r_eval)
        r_sampled_raw[i]  = r_eval

    # -------------------------
    # NUMERICAL STABILITY
    # -------------------------
    Bp_sampled_seed = np.round(Bp_sampled_seed, 7)
    Bp_sampled_raw  = np.round(Bp_sampled_raw, 7)
    r_sampled_seed  = np.round(r_sampled_seed, 7)
    r_sampled_raw   = np.round(r_sampled_raw, 7)

    # -------------------------
    # FILTER VALID ROWS
    # -------------------------
    valid_seed_idx = [
        i for i in range(nrows_seed)
        if not np.allclose(Bp_sampled_seed[i], 0, atol=1e-12)
    ]

    valid_raw_idx = [
        j for j in range(nrows_raw)
        if not np.allclose(Bp_sampled_raw[j], 0, atol=1e-12)
    ]

    print("Valid seed rows:", len(valid_seed_idx))
    print("Valid raw rows:", len(valid_raw_idx))

    if len(valid_seed_idx) == 0 or len(valid_raw_idx) == 0:
        print("No valid rows to match. Skipping galaxy.")
        continue

    # -------------------------
    # BUILD COST MATRIX
    # -------------------------
    cost_matrix = np.zeros((len(valid_seed_idx), len(valid_raw_idx)))

    for ii, i in enumerate(valid_seed_idx):

        seed_bp = Bp_sampled_seed[i]
        seed_r  = r_sampled_seed[i]

        seed_bp_norm = seed_bp / (np.linalg.norm(seed_bp) + 1e-12)

        for jj, j in enumerate(valid_raw_idx):

            raw_bp = Bp_sampled_raw[j]
            raw_r  = r_sampled_raw[j]

            raw_bp_norm = raw_bp / (np.linalg.norm(raw_bp) + 1e-12)

            bp_error = np.sqrt(np.mean((seed_bp_norm - raw_bp_norm) ** 2))
            r_error  = np.sqrt(np.mean((seed_r - raw_r) ** 2))

            cost_matrix[ii, jj] = bp_error + 0.3 * r_error

    # -------------------------
    # GLOBAL MATCHING
    # -------------------------
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    matches = {}

    for ii, jj in zip(row_ind, col_ind):
        seed_i = valid_seed_idx[ii]
        raw_j  = valid_raw_idx[jj]
        matches[seed_i] = raw_j

    print("Matches found:", len(matches))

    # -------------------------
    # SAVE RESULTS
    # -------------------------
    os.makedirs(f"{base}/results/files_corrected", exist_ok=True)

    np.savetxt(
        f"{base}/results/files_corrected/seed_at_diskfractions_{gal_index}.dat",
        Bp_sampled_seed
    )

    np.savetxt(
        f"{base}/results/files_corrected/raw_at_diskfractions_{gal_index}.dat",
        Bp_sampled_raw
    )

    np.savetxt(
        f"{base}/results/files_corrected/seed_row_matches_{gal_index}.txt",
        np.array(list(matches.items())),
        fmt="%d",
        header="SEED_row RAW_row"
    )

    print("Finished galaxy:", gal_index)

print("\nAll galaxies processed.")
