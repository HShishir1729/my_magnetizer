!# Copyright (C) 2018,2019 Luiz Felippe S. Rodrigues, Luke Chamandy
!#
!# This file is part of Magnetizer.
!#
!# Magnetizer is free software: you can redistribute it and/or modify
!# it under the terms of the GNU General Public License as published by
!# the Free Software Foundation, either version 3 of the License, or
!# (at your option) any later version.
!#
!# Magnetizer is distributed in the hope that it will be useful,
!# but WITHOUT ANY WARRANTY; without even the implied warranty of
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!# GNU General Public License for more details.
!#
!# You should have received a copy of the GNU General Public License
!# along with Magnetizer.  If not, see <http://www.gnu.org/licenses/>.
!#
! Global (not galaxy specific) input parameters
! Note to developers: this file should be kept strongly commented, as it will
! will
! also serve as documentation for the parameters
module global_input_parameters
  ! Contains several switches to control the behaviour of the code
  !use messages, only: error_message !CJ
  !use random, only: random_normal   !CJ 
  implicit none

  !lfsr: The following should be elsewhere
  integer :: ngals = -1
  integer :: number_of_redshifts = -1
  integer :: iread=0 !global variable for storing the index of redshift

  ! -------------------------------------------------------
  ! Input and output file settings
  ! -------------------------------------------------------
  ! Separate files for input and output
  logical :: p_IO_separate_output = .true.
  ! If p_IO_separate_output==True, the use the following for the outputfile
  character (len=300) :: output_file_name = 'magnetized_galaxies_output.hdf5'
  character (len=300) :: input_file_name  = 'magnetized_galaxies_input.hdf5'
  ! Chunking options
  logical :: p_IO_chunking = .true.
  integer :: p_IO_number_of_galaxies_in_chunks = 1000
  ! Compression options
  logical :: p_IO_compression = .true. ! requires chunking!
  integer :: p_IO_compression_level = 6
  ! Number of galaxies after which output file is updated
  integer :: p_ncheckpoint = 5000
  ! The master/root process takes part in the calculation? If .true. it will
  ! distribute work and do some work by himself once every
  ! (nproc+int(nproc/p_master_skip)) galaxies.
  ! NB Currently (March/2021), there is a bug in p_master_works_too=.true.
  ! See issue #7 on GitHub for details
  logical :: p_master_works_too = .false.
  double precision :: p_master_skip = 2
  ! Sets which quantities should be added to the output see doc files for
  ! possible values
  character(15), dimension(35) :: output_quantities_list = ' '
  ! Alternatively, tries to output everything.
  logical :: output_everything = .false.

  ! -------------------------------------------------------
  ! Extra debug / custom output switches (HS: seed debugging switch added (17/03/26))
  ! -------------------------------------------------------
  logical :: write_bfloor = .false.
  logical :: write_afloor = .false.
  logical :: write_bp     = .false.
  logical :: write_r      = .false.
  logical :: write_seed   = .false.
  logical :: write_seed_r = .false.
  character(len=100) :: output_suffix = '_fid_100times_random'

  namelist /io_parameters/ &
    output_file_name, input_file_name, p_IO_separate_output, &
    p_IO_chunking, p_IO_number_of_galaxies_in_chunks, p_IO_compression, &
    p_IO_compression_level, &
    p_master_works_too, p_master_skip, p_ncheckpoint, &
    output_quantities_list, output_everything, &
    write_bfloor, write_afloor, write_bp, write_r, & ! HS: Change made for the seed debugging mode (17/03/26)
    write_seed, write_seed_r, &
    output_suffix

  ! -------------------------------------------------------
  ! Run settings and timestepping parameters
  ! -------------------------------------------------------
  character (len=50) :: model_name = 'Magnetized SAM'
  integer :: info = 1
  integer :: nsteps_0 = 1000
  double precision :: p_courant_v = 0.09 !0.4 for Damp=T and nxphys=151
  double precision :: p_courant_eta = 0.09 !0.3 for Damp=T and nxphys=151
  logical :: p_variable_timesteps = .true.
  integer :: p_nsteps_max = 20000
  integer :: p_MAX_FAILS = 6 ! Number of time-step reductions before giving up
  ! Debug mode: all timesteps are included in the output, but
  ! only 1 snapshot is used.
  logical :: p_oneSnaphotDebugMode = .false.
  ! Runs without solving the dynamo equations
  logical :: p_no_magnetic_fields_test_run = .false.
  ! Seed for the random number generator (can be any integer different from 0)
  integer :: p_random_seed = 17
  ! Maximum walltime for the run in seconds (if negative, this is unlimited)
  double precision :: p_max_walltime = -1

  namelist /run_parameters/ &
    model_name, info, nsteps_0, p_courant_v, p_courant_eta, &
    p_variable_timesteps, p_nsteps_max, p_oneSnaphotDebugMode, &
    p_no_magnetic_fields_test_run, p_MAX_FAILS, p_random_seed, &
    p_max_walltime

  ! -------------------------------------------------------
  ! Grid settings
  ! -------------------------------------------------------
  ! If rdisk(t_i) < rdisk(t_{i-1}), rescale the B and alpha WITH the disk
  ! (this keeps the number of grid points and is equivalent to a change of
  !  units)
  logical :: p_rescale_field_for_shrinking_disks=.true.
  ! If rdisk(t_i) > rdisk(t_{i-1}), rescale the B and alpha WITH the disk
  ! (this keeps the number of grid points and is equivalent to a change of
  !  units)
  logical :: p_rescale_field_for_expanding_disks=.false.
  logical :: p_scale_back_f_array = .true.
  ! Reference grid size (used both initially and for storage)
  integer :: p_nx_ref=101  !151
  ! Maximum possible grid size, in the case of varying number of grid points
  ! (this is mostly for debugging)
  integer :: p_nx_MAX=10000
  ! Uses a fixed dimensional grid, with r_max_kpc set using the the maximum
  ! half mass radius the galaxy reaches over the entire History
  logical :: p_use_fixed_physical_grid=.false.
  ! The maximum radius to use for computations divided by the half mass radius
  ! i.e. rmax = p_rmax_over_rdisk * rdisk
  double precision :: p_rmax_over_rdisk = 2.75d0
  ! Which method to use for the interpolation, valid options are:
  ! 'simple' (linear, without using fgsl), 'linear' (fgsl),
  ! 'cubic_spline' (fgsl) and 'akima' (fgsl)
  character(len=100) :: p_interp_method = 'simple'
  !!  Warning: the fgsl interpolators have been temporarily disabled

  namelist /grid_parameters/ &
    p_rescale_field_for_shrinking_disks, p_rescale_field_for_expanding_disks, &
    p_scale_back_f_array, p_nx_ref, p_nx_MAX, p_use_fixed_physical_grid, &
    p_rmax_over_rdisk, p_interp_method

  ! -------------------------------------------------------
  ! Dynamo equations parameters and switches
  ! -------------------------------------------------------
  ! ALPHA QUENCHING
  ! Works with alg_quench=F; Set to T for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Dyn_quench = .true.
  ! Works with dyn_quench=F; Set to T for algebraic alpha quenching; F for no quenching
  logical :: Alg_quench = .false.

  ! CLOSURE APPROXIMATION
  ! Set to T for FOSA, F for minimal tau approximation
  logical :: Damp = .false.

  ! SEED MAGNETIC FIELD
  !Seed field amplitude as a fraction of equipartition magnetic field strength
  double precision :: frac_seed = 0.01d0 !CJ frac_seed = 0.01d0
  ! Selects seed type: 'random', 'decaying' or 'fraction'
  ! fraction -> The seed field is a fixed fraction of the equipartition field
  !             Bseed = Bs = frac_seed*Beq
  ! random   -> The value is drew from a Gaussian distribution with variance Bs
  !             Bseed = Bs*random_gaussian
  ! decaying -> The value decays exponentially with radius
  !             Bseed = r/rmax*(1.d0-r/rmax)**p_nn_seed*dexp(-r/p_r_seed_decay)
  character(len=20) :: p_seed_choice = 'random'
  double precision :: p_r_seed_decay = 15.0d0
  integer:: p_nn_seed = 2 !Only relevant if Rand_seed=F
  ! CEILING ON ALPHA EFFECT
  ! Set to T to put a ceiling for alpha at alpceil*v
  logical :: Alp_ceiling = .true.
  ! ALPHA^2 EFFECT
  ! Set to T to include alpha^2 effect; set to F to use alpha-omega approximation equations
  logical :: Alp_squared= .false.
  ! KRAUSE'S LAW
  ! Set to T for alpha effect to decrease in proportion to omega (Krause's formula)
  logical :: Krause= .true.
  ! ADVECTION
  ! Set to F to turn off advection of the magnetic field
  logical :: Advect= .true.
  ! TURBULENT DIFFUSION
  !Set to F to turn off turbulent diffusion
  logical :: Turb_dif= .true.
  ! Imposes Neumann boundary condition at the maximum radius
  ! I.e. at R=Rmax, symmetric boundary
  logical :: p_neumann_boundary_condition_rmax = .false.
  ! ALPHA EFFECT
  !Factor determining strength of alpha effect
  double precision :: C_alp = 1.0d0  !CJ edit. C_alp=1.0d0
  ! FLOOR
  ! Uses a source term to achieve an effective floor for the magnetic field
  logical :: lFloor = .true.
  !Factor determining strength of the floor
  double precision :: C_floor = 1.0d0
  ! Floor localization (Delta r = p_floor_kappa * p_ISM_turbulent_length)
  ! NB: if using p_use_fixed_turbulent_to_scaleheight_ratio,
  !     this gets inconsistent!
  double precision :: p_floor_kappa = 10d0
  ! Stochastically changes the sign of the floor every Delta r
  logical :: p_space_varying_floor = .true.
  ! Stochastically changes the sign of the floor every
  ! Delta t = p_floor_kappa * tau, where
  !   tau = p_ISM_turbulent_length/(p_ISM_sound_speed_km_s*sqrt(p_ISM_kappa))
  logical :: p_time_varying_floor = .true.
  ! RANDOM MAGNETIC FIELD
  !Strength of the rms random magnetic field brms in units of Beq [brms=fmag*Beq]
  double precision :: fmag = 0.9 !SG: 1.0 in J24, 0.5 in R19 
  ! DIFFUSIVE MAGNETIC HELICITY FLUX
  !Ratio kappa_t/eta_t of turbulent diffusivities of alpha_m and B
  double precision :: R_kappa = 1.5 !SG: 1.5 in J24, 0.3 in R19
  ! in J24 R_kappa = 1.5 SG
  ! in R19 R_kappa = 1.0d0 CJ
  
  double precision :: fBeq =1.0!SG: Default value should be 1.0. (old comments:added by CJ. Default value 0.5)
        ! Used in calculating equipartition field strength Beq = fBeq*sqrt(4*pi*rho)*v (used in profile.f90). 
        ! defined by Eq. 11 of R19. Used by profiles.f90

  ! Fraction of magnetic pressure contribution in the midplane gas pressure.
  ! It can also change the cosmic ray pressure as they are in equipartition.
  double precision:: frac_mag_P = 1.0 ! added by SG! Default value should be 1.0.


  namelist /dynamo_parameters/ &
    Dyn_quench, Alg_quench, lFloor, Damp, &
    frac_seed, p_seed_choice, p_r_seed_decay, p_nn_seed, &
    Alp_ceiling, Alp_squared, Krause, Advect, Turb_dif, &
    p_neumann_boundary_condition_rmax, C_alp, C_floor, &
    p_floor_kappa, p_space_varying_floor, p_time_varying_floor, &
    fmag, R_kappa, fBeq, frac_mag_P


  ! -------------------------------------------------------
  ! Interstellar medium
  ! -------------------------------------------------------
  ! Sound speed (in km/s)
  double precision :: p_ISM_sound_speed_km_s = 25d0 !#CJ p_ISM_sound_speed_km_s= 10d0
  ! Ratio between turbulent velocity and sound speed
  double precision :: p_ISM_kappa = 1d0
  ! Ratio between turbulent pressure and turbulent magnetic field pressure
  ! NB if simplified_pressure is on, this correspond to the ratio between
  ! turbulent pressure and _total_ magnetic field pressure
  double precision :: p_ISM_xi =  0.25!1d0 !not used in current version 
                !due to corrections in pressure_equilibrium.f90 CJ
  ! Adiabatic index of the ISM
  double precision :: p_ISM_gamma = 5d0/3d0
  ! Turbulent length (in kpc)
  double precision :: p_ISM_turbulent_length = 0.1
  ! Uses the scale height as an upper limit for the turbulent scale
  logical :: p_limit_turbulent_scale = .true.

  ! Ratio between stellar scale-height and stellar scale-radius
  double precision :: p_stellarHeightToRadiusScale = 1d0/7.3

  ! Ratio between molecular scale-height and molecular scale-radius
  double precision :: p_molecularHeightToRadiusScale = 0.032

  ! Ratio between total gas scale length and stellar scale length
  double precision :: p_gasScaleRadiusToStellarScaleRadius_ratio = 1d0

  ! The fraction of molecular gas. Default value 1 as in R19 and J24 # SG: added for scaling relation project
  double precision :: frac_mol = 0.0d0!added by SG! 0.5! default value 0.0 in G25

  ! Molecular fraction calculation
  ! Blitz&Rosolowsky alpha
  double precision :: p_Rmol_alpha = 0.92d0
  ! Blitz&Rosolowsky P0 (in erg/cm^3)
  double precision :: p_Rmol_P0 = 4.787d-12

  ! Check whether the hydrostatic equilibrium solution is correct
  logical :: p_check_hydro_solution = .false.

  ! When set to false, this substitutes any positive shear by 0
  ! (NB positive shears can only arise from the regularisation)
  logical :: p_allow_positive_shears = .false.

  ! If true: uses a fixed l/h ratio, setting it with
  ! l = p_turbulent_to_scaleheight_ratio * h
  logical :: p_use_fixed_turbulent_to_scaleheight_ratio = .false.
  double precision :: p_turbulent_to_scaleheight_ratio = 0.25
  ! Uses a (very) simplified calculation for the mid-plane pressure
  ! where P_B + P_b = \xi P_{turb}
  ! (alternatively, P_B uses the actual B from the dynamo calculation).
  logical :: p_simplified_pressure = .true.
  ! Includes the rotation curve correction in the calculation of the midplane pressure
  logical :: p_enable_P2 = .false.
  ! Assumes constant scaleheight for r<rreg
  logical :: p_truncates_within_rreg = .false.
  ! Sets the regularisation radius for the rotation curves
  double precision :: p_rreg_to_rdisk = 0.1

  ! Ignore the DM haloes of satellites for the pressure and rotation curve
  ! calculations
  logical :: p_ignore_satellite_DM_haloes = .false.

  ! Minimum density floor (in g/cm^3)
  double precision :: p_minimum_density = 1d-30 ! g cm^-3

  ! Defines what it means to have a negligible disk
  double precision :: p_rdisk_min=0.5 !kpc
  double precision :: Mgas_disk_min=1d4 ! solar masses
  double precision :: rmin_over_rmax=0.001

  logical :: p_use_Pdm = .true.
  logical :: p_use_Pbulge = .true.
  logical :: p_halo_contraction = .true.
  logical :: p_extra_pressure_outputs = .false.
  ! Keeps P2>0 to avoid possible artifacts which the regularisation may introduce
  ! Usually not necessary.
  logical :: p_P2_workaround = .false.

  ! Legacy, probably unnecessary, options
  ! Adopts a \rho \propto sech^2(z/h) profile for gas and stars
  ! NB Defaul: rho \propto exp(-|z|/h)
  logical :: p_sech2_profile = .false.
  logical :: p_Pdm_numerical = .false.
  logical :: p_Pbulge_numerical = .false.
  


  namelist /ISM_and_disk_parameters/ &
    p_ISM_sound_speed_km_s, p_ISM_kappa, p_ISM_xi, p_ISM_gamma, &
    p_ISM_turbulent_length, p_limit_turbulent_scale, &
    p_stellarHeightToRadiusScale, p_molecularHeightToRadiusScale, &
    p_gasScaleRadiusToStellarScaleRadius_ratio, frac_mol, &
    p_Rmol_alpha, p_Rmol_P0, p_check_hydro_solution, &
    p_allow_positive_shears, p_use_fixed_turbulent_to_scaleheight_ratio, &
    p_turbulent_to_scaleheight_ratio, p_simplified_pressure, &
    p_rreg_to_rdisk, p_rdisk_min, Mgas_disk_min, rmin_over_rmax, &
    p_ignore_satellite_DM_haloes, p_enable_P2, p_truncates_within_rreg, &
    p_minimum_density, p_sech2_profile, p_use_Pdm, p_use_Pbulge, &
    p_halo_contraction, p_extra_pressure_outputs, p_Pdm_numerical, &
    p_Pbulge_numerical, p_P2_workaround


  ! -------------------------------------------------------
  ! Ouflows
  ! -------------------------------------------------------
  ! Outflow calculation ('no_outflow'/'vturb'/'superbubble_simple'/'superbubble/wind')
  character(len=21) :: p_outflow_type = 'no_outflow'
  ! Mechanical luminosity associated with the superbubble (in erg/s)
  double precision :: p_outflow_Lsn = 1d38
  ! Ratio between OB associations (superbubbles) rate and supernovae rate
  double precision :: p_outflow_fOB = 0.7
  ! Ratio between supernovae rate and SFR (in \msun^{-1})
  double precision :: p_outflow_etaSN = 9.4d-3
  ! Lifetime of an OB association (in Myr)
  double precision :: p_tOB = 3
  ! Number of supernovae in 1 OB association
  double precision :: p_N_SN1OB = 40
  ! Density of the hot gas (in g/cm^3)
  double precision :: p_outflow_hot_gas_density = 1.7d-27
  ! Inverse of the depletion timescale of molecular gas (H2+He) in units of Gyr-1
  double precision :: p_outflow_nu0 = 0.5
  ! Exponent in Galform's parametrization of the mass loading
  double precision :: p_outflow_alphahot = -3.2
  ! Velocity scale in Galform's parametrization of the mass loading
  double precision :: p_outflow_Vhot = 425 ! eventually, this should be read
                                           ! from the input parameters file
  namelist /outflow_parameters/ &
    p_outflow_type, &
    p_outflow_Lsn, p_outflow_fOB, p_outflow_etaSN, p_tOB, p_N_SN1OB, &
    p_outflow_hot_gas_density, p_outflow_nu0, p_outflow_alphahot, &
    p_outflow_Vhot

  

  ! -------------------------------------------------------
  ! Observables
  ! -------------------------------------------------------
  ! CR spectral index
  double precision :: p_obs_CR_alpha = 3d0
  ! Dust emissivity power
  ! emissivity is assumed to be proportional to n_H2^alpha
  double precision :: p_obs_dust_alpha = 3d0
  logical :: p_obs_scale_with_z = .true.
  logical :: p_obs_ignore_small_scale = .true.
  ! If set to true, will compute observables for all possible redshifts
  logical :: p_obs_use_all_redshifts = .true.
  ! If p_obs_use_all_redshifts==F, the indices of specific resdhifts
  ! can be chosen using e.g. p_obs_redshift_indices = 45, 44, 43
  integer, dimension(50) :: p_obs_redshift_indices = -1
  ! -------------------------------------------------------
  ! SFR - Star Formation Rate - Currently added to namelist Observables_parameters 
  ! -------------------------------------------------------
  !This is used to correct for SFR offset between galform model and observation
  !Relevant in input_parameters.f90
  integer :: correct_SFR_offset = 0  !1  !CJ edit
 ! Offset factor of ratio between observed SFR and galform SFR.   
 ! The first element corresponds to offset at the highest redshift
 ! Initially this factor is set to 1; will be read from a file if SFR_offset_correction = 1  
  double precision, dimension(50) :: SFR_offset_fac = 1.0 
  character (len=80) :: SFR_offset_file_name = &
      !'/scratch/ccharles/galform_output/Lacey14_SFRcor_offset_B13.txt' !B13 after correcting for SFR
      !'/home/sukanta/galform_output/Lacey14_SFRcor_offset_HB06.txt' !HB06 after correcting for SFR (used in J24)
      !'/home/sukanta/galform_output/Lacey14_SFR_z_offset_MD14.txt' ! Madau & Dickinson 2014
      '/home/sukanta/galform_output/Lacey14_SFR_z_offset_T25.txt'  ! Trina+26 
!      '/scratch/ccharles/galform_output/Lacey14_SFR_offset_B13.txt' !B13
!      '/scratch/ccharles/galform_output/Lacey14_SFR_offset_HB06.txt' !HB06
  double precision :: SFR_th =3.0 !0.1 !added by CJ. Used in CJ et. al. 2023 RLF paper 
                                   !Above SFR_th turbulent speed rises sharply
                                   !Below SFR_th turbulent speed is a constant 
                                   
  integer :: v0_depends_z = 0  !CJ edit. Also need to choose proper c_s model
  double precision, dimension(50) :: v0_z_ay = 1.0 
  character (len=80) :: v0_z_file_name = &
      '/home/sukanta/galform_output/v0_z_50p.txt' !Mean of z-v0 relationfrom Jimnez 2022
  integer :: include_random_B =1 
                ! 0 - Only large scale field
                ! 1 - Isotropic small scale field
                ! 2 - Anisotropic small scale field
  namelist /observables_parameters/ &
    p_obs_CR_alpha, p_obs_dust_alpha, &
    p_obs_ignore_small_scale, p_obs_scale_with_z, &
    p_obs_redshift_indices, p_obs_use_all_redshifts, & 
    correct_SFR_offset, SFR_offset_fac, &
    SFR_offset_file_name, include_random_B, &
    v0_depends_z, v0_z_ay, v0_z_file_name 
     

  contains

  subroutine read_global_parameters(global_pars_filename)
    ! Reads the multiple namelists in global parameters file
    implicit none
    character(len=*), intent(in) :: global_pars_filename
    integer :: u

    open(newunit=u,file=global_pars_filename) ! Fortran 2008.. requires new gfortran
!     u = 17; open(unit=u,file=global_pars_filename) ! Old Fortran
    ! Reads all the namelists
    ! Note: Rewinding makes the order of the namelists in the file unimportant
    read(u,nml=run_parameters); rewind(u)
    read(u,nml=io_parameters); rewind(u)
    read(u,nml=grid_parameters); rewind(u)
    read(u,nml=dynamo_parameters); rewind(u)
    read(u,nml=outflow_parameters); rewind(u)
    read(u,nml=ISM_and_disk_parameters); rewind(u)
    read(u,nml=observables_parameters)
    close(u)
  end subroutine read_global_parameters

!Probably not the best position to put this file


! function to compute ISM turbulent speed. L18 assumes ISM turbulent speed to be
! be a constant (10 km/s). Here it is a function of SFR, Mstar etc.
  function ISM_sound_speed_km_s(Mstar, SFR, Mgas) result(cs)
    double precision, intent(in) :: Mstar, SFR, Mgas
    double precision :: cs, gas_frac, lSFR, sigma_log_cs, log_cs 
    integer :: cs_model

    cs_model = 24 ! 18
              ! 0 is constant cs 
              ! 1 is K16F
              ! 2 is K16G1
              ! 3 is K16G2
              ! 4 is empiracal model 1 - power law fit  to 50th quantile of sfr vs cs data 
              ! 5 is empiracal model 2 - power law fit  to 75th quantile of sfr vs cs data 
              ! 6 is empiracal model 3 - power law fit  to 25th quantile of sfr vs cs data 
              ! 7 is empiracal model 4 - if SFR < 1 Msun, cs= constant. Else cs SFR^alpha 
              ! 8 is stochastic model (stm)  for cs. A Guassian distributed cs around median 
              ! 9 is K16G3. Similar to K16G1.  But cs propto SFR**0.8 
              ! 10 is K16G4. Similar to K16G1.  But cs propto SFR**0.7 
              ! 11 is K16G5. Similar to K16G1.  But cs propto SFR**0.6 
              ! 12 is K16G6. Similar to K16G5.  But cs = 20 at SFR = 3 
              ! 13 is K16G7. Similar to K16G5.  But cs = 20 at SFR = 3. Also cs_max=100 
              ! 14 is K16G8. Similar to K16G5.  But cs = 25 at SFR = 3. After cs propto SFR**0.5 
              ! 15 is K16G9. Similar to K16G5.  But cs = 40 at SFR = 3. After cs propto SFR**0.52 
              ! 16 is K16G10. Similar to K16G5.  But cs = 20 at SFR = 3. After cs propto SFR**0.45 
              ! 17 is J22p50. turbulent speed at different z from Jimnez et al  2022 



!!!!!
      if (cs_model == 0) then ! Constant profile
        cs  = 25.0!10.0 + 2.8*(SFR/SFR_th)**1.0  !10.0
!!!!! 
      else if (cs_model == 1) then 
        cs  = 10.0*SQRT(SFR/1.0)
!!!!! 
      else if (cs_model == 2) then
        cs  = 10.0*(SFR/1.0)
!!!!! 
      else if (cs_model == 3) then
        gas_frac = Mgas/(Mgas+Mstar)
        cs       = 10.0*(SFR/1.0)*(1.0/gas_frac**2.0) 
!!!!! 
      else if (cs_model == 4) then
        lSFR = LOG10(SFR)
        cs   = 10.0**( 0.178*lSFR + 1.56)
!!!!! 
      else if (cs_model == 5) then
        lSFR = LOG10(SFR)
        cs   =  10.0**( 0.17*lSFR + 1.68)
!!!!! 
      else if (cs_model == 6) then
        lSFR = LOG10(SFR)
        cs   = 10.0**( 0.156*lSFR + 1.45)
!!!!! 
      else if (cs_model == 7) then
        lSFR = LOG10(SFR)
        if (lSFR < 0.0) then
           cs   = SQRT(3.0) * 10.0** 1.153
        else
           cs   = SQRT(3.0) * 10.0**( 1.153 + lSFR * 0.3 )
        end if  
!!!!! 
      else if (cs_model == 8) then
        lSFR   = LOG10(SFR)
        log_cs = 0.17*lSFR + 1.3224
       ! cs     = SQRT(3.0) * 10**(0.18*random_normal() + log_cs)
        cs = 10
!!!!! 
      else if (cs_model == 9) then
        cs  = 10.0*(SFR/1.0)**0.8
!!!!! 
      else if (cs_model == 10) then
        cs  = 10.0*(SFR/1.0)**0.7
!!!!! 
      else if (cs_model == 11) then
        cs  = 3.28060192 * (1.0 + (SFR / (6.290e-6)) ** 0.200712215)! Feducial model for G24
!!!!! 
      else if (cs_model == 12) then
        cs  = 20.0*(SFR/3.0)**0.6
        if (cs < 20.0) then
           cs = 20.0
        end if 
!!!!! 
      else if (cs_model == 13) then
        cs  = 20.0*(SFR/3.0)**0.6
        if (cs > 100.0) then
           cs = 100.0
        end if 
        if (cs < 20.0) then
           cs = 20.0
        end if 
!!!!! 
      else if (cs_model == 14) then
        cs  = 27.0*(SFR/4.0)**0.55
        if (cs < 27.0) then
           cs = 27.0
        end if 

!!!!! 
      else if (cs_model == 15) then
        cs  = 41.0*(SFR/4.0)**0.56
        if (cs < 41.0) then
           cs = 41.0
        end if 

!!!!!!!! 
      else if (cs_model == 16) then
        cs  = 18.0*(SFR/4.0)**0.49
        if (cs < 18.0) then
           cs = 18.0
        end if 

!!!!!!!! 
      else if (cs_model == 17) then
        cs  = v0_z_ay(iread) 

!!!!! 
      else if (cs_model == 18) then !J24: Fiducial model (fitted median)
        if (SFR >= SFR_th) then 
           cs  = 20.0*(SFR/SFR_th)**0.25
        else 
           cs = 20.0
        end if 

!!!!!

!!!!! 
      else if (cs_model == 19) then
        if (SFR >= SFR_th) then 
          cs  = 40.0*(SFR/SFR_th)**0.55
        else 
           cs = 40.0
        end if 

!!!!!!!! 
      else if (cs_model == 20) then
        if (SFR >= SFR_th) then
           cs  = 15.0*(SFR/SFR_th)**0.45
        else 
           cs = 15.0
        end if 

!!!!!!!! 
      else if (cs_model == 21) then
        cs  = 25.0
!!!!!!!!      
      else if (cs_model == 22) then
         if (SFR <= 1.0) then 
            cs = 10.0
         else 
            cs = 10.0 * SFR**0.5
            if (cs>= 60) then
               cs = 60.0 
            end if
         end if
!!!!!!!! 

      else if (cs_model == 23) then
        if (SFR > 0.1) then
           cs = 1.256 * (SFR / 9.6502d-3)**0.2416
        else
           cs = 22.0
        end if

     ! else if (cs_model == 23) then
      !  if (SFR >=0.1) then
      !     cs  = 1.256 * (SFR/(9.6502 * 10**(-3)))**0.2416 ! SG
      !  else 
      !     cs = 22.0
      !  end if 

!!!!!!!!
      else if (cs_model == 24) then !SG added!
         if (SFR <= SFR_th) then 
            cs = 25.0
         else 
            cs = (25.0*(SFR/SFR_th)**0.50)
          ! Set the maximum value of cs to 60
!            if (cs > 100.0) then
!               cs = 100.0
!            end if
         end if

!!!!!!!!
      else if (cs_model == 25) then !SG added!
         if (SFR <= SFR_th) then
            cs = 25.0/1.5
         else
           cs = (25.0*(SFR/SFR_th)**0.50)/1.5
          ! Set the maximum value of cs to 60
           if (cs > 60.0) then
               cs = 60.0
           end if
         end if


!!!!!!!!      
      else if (cs_model == 26) then
        if (SFR >= SFR_th) then 
           cs  = 25.0*(SFR/SFR_th)**0.50
        else 
           cs = 25.0
        end if 
        if (cs > 60.0) then
           cs = 60.0
        end if


!!!!!!!! 
      else
        write(*,*) 'set cs_model between given range'
        stop
      end if
!!!!! 
      if (cs < 10.0) then
        cs = 10.0
      end if
!!!!! 

  return
  end function ISM_sound_speed_km_s
    

! Reads the SFR_offset which is ratio between observed SFR and 
! galform SFR at galform output redshifts 
  subroutine read_SFR_offset()
    implicit none
    integer :: u, i

    if (correct_SFR_offset == 1 ) then
      if ( number_of_redshifts  > SIZE(SFR_offset_fac) ) then
         write (*,*) correct_SFR_offset
         write (*,*) number_of_redshifts
         print *,  'Error: SFR and SFR_offset_fac has different size'
         stop
      end if
      open(newunit=u,file=trim(SFR_offset_file_name), status='OLD')
      write (*, *) number_of_redshifts
      do i=1, number_of_redshifts 
         read(u, *) SFR_offset_fac(i)
      end do 
      close(u)
      write (*, *) "SFR offset factor read at all redshifts"
    else
      write (*, *) "SFR offset factor not read. Using galform SFR"
    end if
  end subroutine read_SFR_offset

! Reads the SFR_offset which is ratio between observed SFR and 
! galform SFR at galform output redshifts 
  subroutine read_v0_z()
    implicit none
    integer :: u, i

    if (v0_depends_z == 1 ) then
      if ( number_of_redshifts  > SIZE(v0_z_ay) ) then
         write (*,*) v0_z_ay 
         write (*,*) number_of_redshifts
         print *,  'Error: z and v0 has different size'
         stop
      end if
      open(newunit=u,file=trim(v0_z_file_name), status='OLD')
      write (*, *) number_of_redshifts
      do i=1, number_of_redshifts
         read(u, *) v0_z_ay(i)
      end do
      close(u)
      write (*, *) "v0 read at all redshifts"
    else
      write (*, *) "v0 not a function of z"
    end if
  end subroutine read_v0_z 



end module global_input_parameters
