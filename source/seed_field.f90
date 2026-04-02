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
module seed_field  !Set initial conditions
  use global_input_parameters
  use grid
  use random
  use profiles, only: l, h, Beq, delta_r_floor
  use floor_field
  implicit none

  contains
    subroutine init_seed_field(f)
      integer :: i,var
      double precision, dimension(nx) :: Bseed
      double precision, dimension(:,:), intent(inout) :: f
      double precision :: r1
      ! Initializes the seed field to a fraction of the equipartition field
      Bseed=frac_seed*Beq
      ! Initializes f
      f(:,:)=0.d0

      select case (trim(p_seed_choice))
        case('fraction')
          ! The field is a fixed fraction of the equipartition field
          f(:,1)=-Bseed
          f(:,2)= Bseed
        case('decaying')
          ! The field
          r1 = p_r_seed_decay/r_max_kpc
          f(:,1)=-Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
          f(:,2)= Bseed*r*(1.d0-r)**p_nn_seed*dexp(-r/r1)
        case('random')
          ! The field at each point (and component) is a Gaussian drawing from
          ! a distribution with variance Bseed
          do var=1,2 !Seed r and phi components of B
            do i=1,nx
              f(i,var)=Bseed(i)*random_normal()
            enddo
          enddo
        case('floor')
          ! Uses floor target field as seed (leave Br=0)
          Bseed = frac_seed*compute_floor_target_field(r,l,h,Beq,delta_r_floor)
          f(:,1)= 0d0
          f(:,2)= Bfloor_sign*Bseed
        case default
          print *, 'init_seed: Invalid option ',p_seed_choice
          stop
      end select

      call write_seed_to_file(f)

    end subroutine init_seed_field

    subroutine write_seed_to_file(f)

      use input_params
      use global_input_parameters
      use grid, only: r

      implicit none
      double precision, intent(in) :: f(:,:)
      integer :: i, unit, unit_r
      character(len=300) :: seed_fname, r_fname
    
!      open(newunit=unit, file='/home/shishir/test/B_op_fiducial_single/seed_output_6147fid.dat', &
!           status='unknown', position='append', action='write')
!      write(unit,'(*(E20.8,1X))') (f(i,2), i=1,size(f,1))
!      close(unit)
      ! File for seed field
      if (write_seed) then
        write(seed_fname,'(A,I0,A,A)') &
             '/home/shishir/test/B_op_fiducial_single/seed/seed_output_', &
             current_gal_id, trim(output_suffix), '.dat'
        open(newunit=unit, file=trim(seed_fname), &
             status='replace', action='write')
        write(unit,'(*(E20.8,1X))') (f(i,2), i=1,size(f,1))
        close(unit)
      endif

      ! File for radius
      if (write_seed_r) then
        write(r_fname,'(A,I0,A,A)') &
             '/home/shishir/test/B_op_fiducial_single/r_seed/radius_seed_', &
             current_gal_id,trim(output_suffix), '.dat'
        open(newunit=unit_r, file=trim(r_fname), &
             status='replace', action='write')
        write(unit_r,'(*(E20.8,1X))') (r(i), i=1,size(r))
        close(unit_r)
      endif

    end subroutine write_seed_to_file
  
end module seed_field
