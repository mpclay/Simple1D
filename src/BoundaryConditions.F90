! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>
!
!> @file BoundaryConditions.F90
!> @author Matthew Clay
!> @brief Boundary conditions for the 1D advection code.
!!
!! Right now this module only implements periodic boundary conditions, but
!! others should be implemented in the future for further investigations.
MODULE BoundaryConditions_m

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   ! Different types of boundary conditions.
   !
   !> Periodic BCs.
   INTEGER(KIND=IWP),PARAMETER,PUBLIC :: PERIODIC = 1_IWP

   ! Module procedures.
   PUBLIC :: BoundaryConditions

CONTAINS

   !> Procedure to set the desired type of boundary conditions.
   !!
   !> @param[in,out] u Solution array.
   !> @param[in] n Number of finite volume cells.
   !> @param[in] ngc Number of ghost cells.
   !> @param[in] bcType Type of boundary conditions to apply.
   SUBROUTINE BoundaryConditions(u, n, ngc, bcType)
      IMPLICIT NONE
      ! Calling arguments.
      REAL(KIND=RWP),DIMENSION(1-ngc:n+ngc),INTENT(INOUT) :: u
      INTEGER(KIND=IWP),INTENT(IN) :: n, ngc, bcType

      ! Set the boundary conditions depending on the desired type.
      SELECT CASE (bcType)
         CASE (PERIODIC)
            u(1-ngc:0) = u(n-ngc+1:n)
            u(n+1:n+ngc) = u(1:ngc)
         CASE DEFAULT
            WRITE(*,*) 'Invalid option for the boundary conditions. Halting.'
            STOP
      END SELECT
   END SUBROUTINE BoundaryConditions

END MODULE BoundaryConditions_m

