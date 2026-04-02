import numpy as np
import h5py
from scipy.interpolate import interp1d
import sys

galaxy_list_file = sys.argv[1]
task_id = int(sys.argv[2])

galaxy_indices = np.loadtxt(galaxy_list_file, dtype=int)
gal_index = galaxy_indices[task_id]

print("Processing galaxy:", gal_index)

base = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new"

raw_bp_file = f"{base}/bp/Bp_profiles_{gal_index}_fid_floor.dat"
raw_r_file  = f"{base}/r/r_profiles_{gal_index}_fid_floor.dat"

h5_file = f"{base}/Lacey14_G25_fiducial_trial_1_corr_single_10galx_floor.hdf5"


def read_variable_rows(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.strip():
                data.append(np.array(line.split(), dtype=float))
    return data


# ---------------------------
# Read raw profiles
# ---------------------------
Bp_profiles = read_variable_rows(raw_bp_file)
r_profiles  = read_variable_rows(raw_r_file)

nrows = len(r_profiles)

Bp_sampled_raw = np.zeros((nrows, 5))
r_sampled_raw  = np.zeros((nrows, 5))

for i in range(nrows):

    r = r_profiles[i]
    Bp = Bp_profiles[i]

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


# ---------------------------
# Read HDF5 data
# ---------------------------
with h5py.File(h5_file, "r") as f:
    Bp = f["Output/Bp"][gal_index - 1]
    r  = f["Output/r"][gal_index - 1]

nr, nz = Bp.shape

Bp_sampled_h5 = np.zeros((nz, 5))
r_sampled_h5  = np.zeros((nz, 5))

for j in range(nz):

    r_profile = r[:, j]
    Bp_profile = Bp[:, j]

    zero_idx = np.where(np.isclose(Bp_profile, 0))[0]

    if len(zero_idx) < 2:
        continue

    start = zero_idx[0]
    end   = zero_idx[-1]

    r_valid = r_profile[start:end+1]
    Bp_valid = Bp_profile[start:end+1]

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

    Bp_sampled_h5[j] = interp_func(r_eval)
    r_sampled_h5[j]  = r_eval


# ---------------------------
# Improved matching algorithm
# ---------------------------
initial_tol = 1e-5
max_attempts = 5

matches = {}

global_best_error = np.inf
global_best_pair = None

for i in range(Bp_sampled_h5.shape[0]):

    row_h5 = Bp_sampled_h5[i]

    if np.allclose(row_h5, 0, atol=1e-12):
        continue

    tol = initial_tol
    best_match = None
    best_error = np.inf

    # Exact matching with tolerance growth
    for attempt in range(max_attempts):

        candidate_matches = []

        for j in range(Bp_sampled_raw.shape[0]):

            row_raw = Bp_sampled_raw[j]

            if row_h5.shape != row_raw.shape:
                continue

            if np.allclose(row_h5, row_raw, atol=tol):
                error = np.sqrt(np.mean((row_h5 - row_raw) ** 2))
                candidate_matches.append((j, error))

        if candidate_matches:
            candidate_matches.sort(key=lambda x: x[1])
            best_match, best_error = candidate_matches[0]
            matches[i] = best_match

            if best_error < global_best_error:
                global_best_error = best_error
                global_best_pair = (i, best_match)

            break

        tol *= 10

    # Fallback matching
    if best_match is None:

        tol = initial_tol

        for attempt in range(max_attempts):

            candidate_matches = []

            for j in range(Bp_sampled_raw.shape[0]):

                row_raw = Bp_sampled_raw[j]

                if np.allclose(row_raw, 0, atol=1e-12):
                    continue

                diff = np.abs(row_h5 - row_raw)
                n_match = np.sum(diff <= tol)

                if n_match >= 3:
                    error = np.sqrt(np.mean((row_h5 - row_raw) ** 2))
                    candidate_matches.append((j, error))

            if candidate_matches:
                candidate_matches.sort(key=lambda x: x[1])
                best_match, best_error = candidate_matches[0]
                matches[i] = best_match

                if best_error < global_best_error:
                    global_best_error = best_error
                    global_best_pair = (i, best_match)

                break

            tol *= 10


# ---------------------------
# Save outputs
# ---------------------------
np.savetxt(
    f"{base}/results/files_corrected/Bp_at_diskfractions_{gal_index}.dat",
    Bp_sampled_raw
)

np.savetxt(
    f"{base}/results/files_corrected/r_at_diskfractions_{gal_index}.dat",
    r_sampled_raw
)

np.savetxt(
    f"{base}/results/files_corrected/Bp_at_diskfractions_{gal_index}_h5.dat",
    Bp_sampled_h5
)

np.savetxt(
    f"{base}/results/files_corrected/r_at_diskfractions_{gal_index}_h5.dat",
    r_sampled_h5
)

np.savetxt(
    f"{base}/results/files_corrected/row_matches_reduced_{gal_index}.txt",
    np.array(list(matches.items())),
    fmt="%d",
    header="H5_row RAW_row"
)

if global_best_pair is not None:
    print("Best overall match:", global_best_pair, "Error:", global_best_error)

print("Finished galaxy:", gal_index)
