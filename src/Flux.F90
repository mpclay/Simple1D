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
!> @file Flux.F90
!> @author Matthew Clay
!> @brief Calculate the RHS of the finite volume scheme.
!!
!! Here we use WENO reconstruction to obtain primitives on the cell faces, and 
!! then calculate a numerical flux with the Lax-Friedrichs flux function.
!!
!! In this code we do not do anything fancy to save on computational performance
!! or efficiency. The emphasis is purely on getting WENO correct.
MODULE Flux_m

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   ! Various reconstruction constants.
   !
   !> Reconstruction constants for k = 2 in r, j ordering.
   REAL(KIND=RWP),DIMENSION(3,2),TARGET,PRIVATE :: crjO2 = &
      TRANSPOSE(RESHAPE([1.5_RWP, -0.5_RWP, &
                         0.5_RWP, 0.5_RWP, &
                         -0.5_RWP, 1.5_RWP], [2,3]))
   !> Reconstruction constants for k = 3 in r, j ordering.
   REAL(KIND=RWP),DIMENSION(4,3),TARGET,PRIVATE :: crjO3 = &
      TRANSPOSE(RESHAPE([11.0_RWP/6.0_RWP, -7.0_RWP/6.0_RWP, 1.0_RWP/3.0_RWP, &
                         1.0_RWP/3.0_RWP, 5.0_RWP/6.0_RWP, -1.0_RWP/6.0_RWP, &
                         -1.0_RWP/6.0_RWP, 5.0_RWP/6.0_RWP, 1.0_RWP/3.0_RWP, &
                         1.0_RWP/3.0_RWP, -7.0_RWP/6.0_RWP, 11.0_RWP/6.0_RWP], &
                         [3,4]))
   !> Smooth solution reconstruction constants.
   REAL(KIND=RWP),DIMENSION(:),ALLOCATABLE,PRIVATE :: d
   !> Smoothness indications.
   REAL(KIND=RWP),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: b
   !> Unnormalized nonlinear weights.
   REAL(KIND=RWP),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: a, aH
   !> Normalized nonlinear weights.
   REAL(KIND=RWP),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: w, wH
   !> Epsilon for the nonlinear weights.
   REAL(KIND=RWP),PARAMETER,PRIVATE :: eps = 1.0e-6
   !> Main pointer set to the correct reconstruction constants.
   REAL(KIND=RWP),DIMENSION(:,:),POINTER,PRIVATE :: crj => NULL()

   ! Module procedures.
   PUBLIC :: SetReconstructionOrder, Flux

CONTAINS

   !> Set the various reconstruction constants depending on the order.
   !!
   !> @param[in] n Number of finite volume cells.
   !> @param[in] k Desired polynomial order.
   SUBROUTINE SetReconstructionOrder(n, k)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n, k

      INTEGER(KIND=IWP) :: r, j

      ! Set the reconstruction constants.
      SELECT CASE (k)
         CASE (2)
            ! Set the pointer with the appropriate starting and ending indices
            ! to match what is in the lecture notes.
            crj(-1:1,0:1) => crjO2

            ! Allocate the smooth reconstruction constants.
            ALLOCATE(d(0:1))
            d(0) = 2.0_RWP/3.0_RWP
            d(1) = 1.0_RWP - d(0)

            ! Allocate the smoothness indicators.
            ALLOCATE(b(n,0:1))

            ! Allocate the unnormalized nonlinear weights.
            ALLOCATE(a(n,0:1))
            ALLOCATE(aH(n,0:1))

            ! Allocate the normalized nonlinear weights.
            ALLOCATE(w(n,0:1))
            ALLOCATE(wH(n,0:1))
         CASE (3)
            ! Set the pointer to the appropriate starting and ending indices to
            ! match what is in the lecture notes.
            crj(-1:2,0:2) => crjO3

            ! Allocate the smooth reconstruction constants.
            ALLOCATE(d(0:2))
            d(0) = 3.0_RWP/10.0_RWP
            d(1) = 3.0_RWP/5.0_RWP
            d(2) = 1.0_RWP - d(0) - d(1)

            ! Allocate the smoothness indicators.
            ALLOCATE(b(n,0:2))

            ! Allocate the unnormalized nonlinear weights.
            ALLOCATE(a(n,0:2))
            ALLOCATE(aH(n,0:2))

            ! Allocate the normalized nonlinear weights.
            ALLOCATE(w(n,0:2))
            ALLOCATE(wH(n,0:2))
         CASE DEFAULT
            WRITE(*,*) 'Invalid polynomial order. Halting.'
            STOP
      END SELECT
   END SUBROUTINE SetReconstructionOrder

   !> Routine to calculate the RHS of the finite volume scheme.
   !!
   !! Some clarification about the naming conventions used in the reconstruction
   !! and flux arrays used below:
   !!
   !!    - When looping over each CELL we reconstruct primitives on its right
   !!      and left faces. These are stored in uCellR and uCellL, respectively.
   !!    - When looping over each FACE, we calculate a flux using left and right
   !!      extrapolated values as seen by the FACE. In other words, the flux
   !!      at a face calculated with values extrapolated from the left are
   !!      considered the right-reconstructed values as seen by the cell to the
   !!      left of that face. Vice versa for the right flux, i.e., it is
   !!      calculated with left-reconstructed values as seen by the cell to the
   !!      right of that face.
   !!
   !! TODO:
   !!    1. Get rid of the flux copying at the beginning and end of the domain
   !!       since that enforces periodic boundary conditions.
   !!
   !> @param[in] u Current state vector.
   !> @param[in] wSpeed Wave speed for the scalar equation.
   !> @param[in] n Number of finite volume cells.
   !> @param[in] k Order of the polynomial reconstruction.
   !> @param[in] ngc Number of ghost cells.
   !> @param[in] dx Local grid size.
   !> @param[in,out] Lu RHS of the system.
   SUBROUTINE Flux(u, wSpeed, n, k, ngc, dx, Lu)
      IMPLICIT NONE
      ! Calling arguments.
      REAL(KIND=RWP),DIMENSION(1-ngc:n+ngc),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(IN) :: wSpeed, dx
      INTEGER(KIND=IWP),INTENT(IN) :: n, k, ngc
      REAL(KIND=RWP),DIMENSION(n),INTENT(OUT) :: Lu
      ! Local variables.
      ! Left reconstructed states for each cell using each stencil.
      REAL(KIND=RWP),DIMENSION(1:n,0:k-1) :: uCellL
      ! Right reconstructed states for each cell using each stencil.
      REAL(KIND=RWP),DIMENSION(1:n,0:k-1) :: uCellR
      ! Primitives extrapolated from the left for a given cell face.
      REAL(KIND=RWP),DIMENSION(0:n) :: uL
      ! Primitives extrapolated from the right for a given cell face.
      REAL(KIND=RWP),DIMENSION(0:n) :: uR
      ! Flux calculated with left-extrapolated values.
      REAL(KIND=RWP),DIMENSION(0:n) :: fL
      ! Flux calculated with right-extrapolated values.
      REAL(KIND=RWP),DIMENSION(0:n) :: fR
      ! Upwinded flux.
      REAL(KIND=RWP),DIMENSION(0:n) :: fU
      ! Sums of a and aH for normalization.
      REAL(KIND=RWP) :: aSum, aHSum
      ! Looping index for cells.
      INTEGER(KIND=IWP) :: i
      ! Looping index for stencils.
      INTEGER(KIND=IWP) :: r
      ! Looping index for polynomial reconstruction coefficients.
      INTEGER(KIND=IWP) :: j

      ! Zero out arrays that are incremented.
      uR(:) = 0.0_RWP
      uL(:) = 0.0_RWP
      uCellL(:,:) = 0.0_RWP
      uCellR(:,:) = 0.0_RWP

      ! Figure out the weights based on the polymial order.
      !
      ! 1. Calculate the smoothness indicators, which depend on the stencil size
      !    being used.
      SELECT CASE (k)
         CASE (2)
            DO i = 1, n
               b(i,0) = (u(i+1) - u(i))**2
               b(i,1) = (u(i) - u(i-1))**2
            END DO
         CASE (3)
            DO i = 1, n
               b(i,0) = 13.0_RWP/12.0_RWP*(u(i) - 2.0_RWP*u(i+1) + u(i+2))**2 +&
                        0.25_RWP*(3.0_RWP*u(i) - 4.0_RWP*u(i+1) + u(i+2))**2
               b(i,1) = 13.0_RWP/12.0_RWP*(u(i-1) - 2.0_RWP*u(i) + u(i+1))**2 +&
                        0.25_RWP*(u(i-1) - u(i+1))**2
               b(i,2) = 13.0_RWP/12.0_RWP*(u(i-2) - 2.0_RWP*u(i-1) + u(i))**2 +&
                        0.25_RWP*(u(i-2) - 4.0_RWP*u(i-1) + 3.0_RWP*u(i))**2
            END DO
      END SELECT
      !
      ! 2. Calculate the unnormalized nonlinear weights.
      DO i = 1, n
         DO r = 0, k-1
            a(i,r) = d(r)/(eps + b(i,r))**2
            aH(i,r) = d(k-1-r)/(eps + b(i,r))**2
         END DO
      END DO
      !
      ! 3. Finally, calculate the normalized nonlinear weights.
      DO i = 1, n
         aSum = SUM(a(i,:))
         aHSum = SUM(aH(i,:))
         DO r = 0, k-1
            w(i,r) = a(i,r)/aSum
            wH(i,r) = aH(i,r)/aHSum
         END DO
      END DO
   
      ! For each cell (i), go through each stencil (r) and calculate the k
      ! reconstructed values of u at the left and right faces.
      DO i = 1, n
         DO r = 0, k-1
            DO j = 0, k-1
               ! Obtain the values u_(i+1/2)^r to kth order accuracy.
               uCellR(i,r) = uCellR(i,r) + crj(r,j)*u(i-r+j)

               ! Obtain the values u_(i-1/2)^r to kth order accuracy.
               uCellL(i,r) = uCellL(i,r) + crj(r-1,j)*u(i-r+j)
            END DO
         END DO
      END DO

      ! Now that we have u_(i+-1/2) evaluated using k different stencils,
      ! combine them with the nonlinear weights to obtain an even higher order
      ! reconstruction in smooth regions.
      DO i = 1, n
         DO r = 0, k-1
            uL(i) = uL(i) + w(i,r)*uCellR(i,r)
            uR(i-1) = uR(i-1) + wH(i,r)*uCellL(i,r)
         END DO
      END DO

      ! NOTE: the following flux copies will have to go for non-periodic BCs.
      !
      ! Copy the left-extrapolated values at the right of the domain to the left
      ! face at the beginning of the domain.
      uL(0) = uL(n)
      !
      ! Copy the right-extrapolated values at the left of the domain to the
      ! right face at the end of the domain.
      uR(n) = uR(0)

      ! NOTE: the following flux function has to change for Burgers' equation.
      !
      ! At each face, calculate the fluxes based on the left and right
      ! extrapolated values.
      DO i = 0, n
         fL(i) = wSpeed*uL(i)
         fR(i) = wSpeed*uR(i)
      END DO

      ! NOTE: the wave speed used here will have to change for Burgers' eqn.
      !
      ! Combine the separate fluxes on each face with Lax-Friedrichs.
      DO i = 0, n
         fU(i) = 0.5_RWP*(fL(i) + fR(i) - wSpeed*(uR(i) - uL(i)))
      END DO

      ! For each cell, calculate its increment based on the flux difference. 
      ! Remember, for each cell, a flux array indexed at the cell index is for
      ! the right face, and the flux array indexed at one less than the cell
      ! index is for the left face.
      DO i = 1, n
         Lu(i) = -(fU(i) - fU(i-1))/dx
      END DO
   END SUBROUTINE Flux

   !> Finalizer for the module.
   SUBROUTINE FluxFinalize()
      IMPLICIT NONE
      DEALLOCATE(d)
      DEALLOCATE(b)
      DEALLOCATE(a)
      DEALLOCATE(aH)
      DEALLOCATE(w)
      DEALLOCATE(wH)
      crj => NULL()
   END SUBROUTINE FluxFinalize

END MODULE Flux_m

