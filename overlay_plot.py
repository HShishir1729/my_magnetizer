import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import os

# ================= FILE PATHS =================
file1 = '/home/shishir/test/B_op_fiducial_single_10galx_floor_new/Lacey14_G25_fiducial_trial_1_corr_single_10galx_floor.hdf5'
file2 = '/home/shishir/test/B_op_fiducial_single_10galx_100times_floor_new/Lacey14_G25_fiducial_trial_1_corr_single_10galx_100times_floor.hdf5'

txt_file = '/home/shishir/test/Magnetizer/params_template/galaxy_list10.txt'

# Seed directories
seed_dir_1 = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new/seed"
rseed_dir_1 = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new/r_seed"
beq_dir_1 = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new/seed_beq"

seed_dir_2 = "/home/shishir/test/B_op_fiducial_single_10galx_100times_floor_new/seed"
rseed_dir_2 = "/home/shishir/test/B_op_fiducial_single_10galx_100times_floor_new/r_seed"
beq_dir_2 = "/home/shishir/test/B_op_fiducial_single_10galx_100times_floor_new/seed_beq"

# ================= HELPER =================
def get_profile(arr):
    return arr[0, :] if arr.ndim == 2 else arr

# ================= LOAD GALAXIES =================
galaxies = np.loadtxt(txt_file, dtype=int) - 1

# ================= LOAD HDF5 =================
with h5py.File(file1, 'r') as hdf:
    bp_array_1 = hdf['Output/Br'][()]
    r_array_1 = hdf['Output/r'][()]
    num_z_indices = bp_array_1.shape[2]

with h5py.File(file2, 'r') as hdf:
    bp_array_2 = hdf['Output/Br'][()]
    r_array_2 = hdf['Output/r'][()]

# ================= REDSHIFT (FIXED) =================
redshift_values = np.linspace(6, 0, num_z_indices)

# ================= COLORMAPS =================
cmap1 = plt.cm.winter.reversed()
cmap2 = plt.cm.autumn

norm = Normalize(vmin=0, vmax=6)

sm1 = ScalarMappable(cmap=cmap1, norm=norm)
sm1.set_array([])

sm2 = ScalarMappable(cmap=cmap2, norm=norm)
sm2.set_array([])

# ================= OUTPUT =================
output_pdf = "/home/shishir/test/B_op_fiducial_single_10galx_floor_new/results/code/overlay/galaxy_profiles_overlay_with_seed_zoom_corrected.pdf"

with PdfPages(output_pdf) as pdf:

    for galaxy_index in galaxies:

        fig, (ax, ax_zoom) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

        # ================= MAIN PROFILES =================
        for z_index in range(num_z_indices):

            z_val = redshift_values[z_index]
            color1 = cmap1(norm(z_val))
            color2 = cmap2(norm(z_val))

            # ---------- FILE 1 ----------
            prof1 = bp_array_1[galaxy_index, :, z_index]
            if np.any(~np.isclose(prof1, -99999.003, atol=1)):
                if not (np.all(~np.isfinite(prof1)) or np.nanmax(np.abs(prof1)) > 1e10):
                    r_half = np.max(r_array_1[galaxy_index, :, z_index]) / 2.7
                    x = r_array_1[galaxy_index, :, z_index] / r_half
                    ax.plot(x, prof1, color=color1, linestyle='-', alpha=0.9)
                    ax_zoom.plot(x, prof1, color=color1, linestyle='-', alpha=0.9)

            # ---------- FILE 2 ----------
            prof2 = bp_array_2[galaxy_index, :, z_index]
            if np.any(~np.isclose(prof2, -99999.003, atol=1)):
                if not (np.all(~np.isfinite(prof2)) or np.nanmax(np.abs(prof2)) > 1e10):
                    r_half = np.max(r_array_2[galaxy_index, :, z_index]) / 2.7
                    x = r_array_2[galaxy_index, :, z_index] / r_half
                    ax.plot(x, prof2, color=color2, linestyle=':', alpha=0.9)
                    ax_zoom.plot(x, prof2, color=color2, linestyle=':', alpha=0.9)

        # ================= SEED COLOR (z = 6) =================
        color_z6_file1 = cmap1(norm(6))
        color_z6_file2 = cmap2(norm(6))

        # ================= SEED + BEQ =================
        gal_id = galaxy_index + 1
        seed_vals_all = []

        try:
            # ----- FILE 1 -----
            seed_file_1 = f"{seed_dir_1}/seed_output_{gal_id}_fid_floor.dat"
            r_file_1 = f"{rseed_dir_1}/radius_seed_{gal_id}_fid_floor.dat"
            beq_file_1 = f"{beq_dir_1}/seed_beq_output_{gal_id}_fid_floor.dat"

            if os.path.exists(seed_file_1):
                seed_1 = get_profile(np.loadtxt(seed_file_1))
                r_seed_1 = np.loadtxt(r_file_1)
                beq_1 = get_profile(np.loadtxt(beq_file_1))

                r_norm_1 = r_seed_1 / (np.max(r_seed_1) / 2.7)

                seed_vals_all.extend(seed_1)
#                seed_vals_all.extend(beq_1/100)

                ax.plot(r_norm_1, seed_1, color=color_z6_file1, linestyle='-', linewidth=2)
                ax.plot(r_norm_1, beq_1/100, color='black', linestyle='-', linewidth=2)

                ax_zoom.plot(r_norm_1, seed_1, color=color_z6_file1, linestyle='-', linewidth=2)
                ax_zoom.plot(r_norm_1, beq_1/100, color='black', linestyle='-', linewidth=2)

            # ----- FILE 2 -----
            seed_file_2 = f"{seed_dir_2}/seed_output_{gal_id}_fid_floor.dat"
            r_file_2 = f"{rseed_dir_2}/radius_seed_{gal_id}_fid_floor.dat"
            beq_file_2 = f"{beq_dir_2}/seed_beq_output_{gal_id}_fid_floor.dat"

            if os.path.exists(seed_file_2):
                seed_2 = get_profile(np.loadtxt(seed_file_2))
                r_seed_2 = np.loadtxt(r_file_2)
                beq_2 = get_profile(np.loadtxt(beq_file_2))

                r_norm_2 = r_seed_2 / (np.max(r_seed_2) / 2.7)

                seed_vals_all.extend(seed_2)
#                seed_vals_all.extend(beq_2/100)

                ax.plot(r_norm_2, seed_2, color=color_z6_file2, linestyle=':', linewidth=2)
                ax.plot(r_norm_2, beq_2/100, color='black', linestyle=':', linewidth=2)

                ax_zoom.plot(r_norm_2, seed_2, color=color_z6_file2, linestyle=':', linewidth=2)
                ax_zoom.plot(r_norm_2, beq_2/100, color='black', linestyle=':', linewidth=2)

        except Exception as e:
            print(f"Skipping seed/beq for galaxy {gal_id}: {e}")

        # ================= ZOOM SCALING =================
        if len(seed_vals_all) > 0:
            vals = np.array(seed_vals_all)
            vals = vals[np.isfinite(vals)]
            if len(vals) > 0:
                ymax = np.max(np.abs(vals)) * 2
                ax_zoom.set_ylim(-ymax, ymax)

        # ================= LABELS =================
        ax.set_ylabel(r"$B_\phi\ (\mu{\rm G})$", fontsize=12)
        ax_zoom.set_ylabel("Zoom (Seed scale)", fontsize=11)
        ax_zoom.set_xlabel(r"$r / (r_{1/2})$", fontsize=12)

        ax.set_title(f"Galaxy Index: {gal_id}", fontsize=13)

        # ================= COLORBARS =================
        cbar1 = fig.colorbar(sm1, ax=ax, fraction=0.046, pad=0.04)
        cbar1.set_label("Redshift (Random)", fontsize=10)
        cbar1.ax.invert_yaxis()

        cbar2 = fig.colorbar(sm2, ax=ax, fraction=0.046, pad=0.10)
        cbar2.set_label("Redshift (100× Random)", fontsize=10)
        cbar2.ax.invert_yaxis()

        # ================= LEGEND =================
        legend_lines = [
            Line2D([0], [0], color='black', linestyle='-', label='Random'),
            Line2D([0], [0], color='black', linestyle=':', label='100× Random'),
            Line2D([0], [0], color='black', linestyle='-', linewidth=2, label='Beq (Random)'),
            Line2D([0], [0], color='black', linestyle=':', linewidth=2, label='Beq (100× Random)')
        ]
        ax.legend(handles=legend_lines, loc='upper right')

        ax.grid(alpha=0.3)
        ax_zoom.grid(alpha=0.3)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

print(f"Saved multipage PDF: {output_pdf}")
