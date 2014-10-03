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
!> @file Simple1D.F90
!> @author Matthew Clay
!> @brief 1D finite volume code for solving advection and Burgers' equations.
!!
!! Following the lecture notes of Dr. Shu ("Essentially Non-oscillatory and
!! Weighted Essentially Non-oscillatory Schemes for Hyperbolic Conservation
!! Laws," Chapter 4 in "Lecture Notes in Mathematics" 1998), we write a 1D
!! finite volume code to solve the 1D scalar advection and Burgers' equations.
!! Primitives are reconstructed at the cell faces with WENO reconstruction, and
!! a Lax-Friedrichs flux is used to evaluate the numerical flux. The schemes
!! currently available include:
!!
!! Time:
!!
!!    1. 3rd order TVD RK scheme of Shu.
!!
!! Space:
!!
!!    1. 3rd order FV WENO.
!!    2. 5th order FV WENO.
!!
!! TODO:
!!
!!    1. Need to put the user inputs into an input file.
PROGRAM Simple1D

   ! Required modules.
   USE Parameters_m,ONLY: RWP, IWP, PI
   USE BoundaryConditions_m,ONLY: BoundaryConditions, PERIODIC
   USE Flux_m,ONLY: SetReconstructionOrder, Flux, FluxFinalize

   IMPLICIT NONE

   ! The different systems that can be solved with this code. Note that this
   ! point only ADVECTION is implemented.
   !
   !> Advection equation.
   INTEGER(KIND=IWP),PARAMETER :: ADVECTION = 0_IWP
   !> Burgers' equation.
   INTEGER(KIND=IWP),PARAMETER :: BURGERS = 1_IWP
   !
   ! The different initial conditions that can be used.
   !
   !> Sinusoidal.
   INTEGER(KIND=IWP),PARAMETER :: SINUSOIDAL = 0_IWP
   !> Square wave.
   INTEGER(KIND=IWP),PARAMETER :: SQUARE_WAVE = 1_IWP

   ! User inputs for the simulation.
   !
   !> Governing equation to solve.
   !! NOTE: Currently only scalar advection is implemented.
   INTEGER(KIND=IWP),PARAMETER :: sys = ADVECTION
   !> Length of the domain.
   REAL(KIND=RWP),PARAMETER :: L = 2.0_RWP
   !> Start of the domain.
   REAL(KIND=RWP),PARAMETER :: xstart = 0.0_RWP
   !> Number of finite volume cells.
   INTEGER(KIND=IWP),PARAMETER :: n = 32_IWP
   !> Order of polynomial reconstruction. Options are 2 or 3.
   INTEGER(KIND=IWP),PARAMETER :: k = 3_IWP
   !> CFL for time integration. Used when dynamically evaluating time step.
   REAL(KIND=RWP),PARAMETER :: cfl = 0.01_RWP
   !> Time step for integration. Only used when this is positive.
   REAL(KIND=RWP),PARAMETER :: dtConst = 0.000001_RWP
   !> End time for the simulation. Used when nEnd is negative.
   REAL(KIND=RWP),PARAMETER :: tend = 2.0_RWP
   !> Total number of iterations to perform. Only used when this is positive.
   INTEGER(KIND=IWP),PARAMETER :: nEnd = 2000000_IWP
   !> Type of initial conditions.
   INTEGER(KIND=IWP),PARAMETER :: ics = SINUSOIDAL
   !> Wave speed for the system if the scalar advection equation is being used.
   REAL(KIND=RWP),PARAMETER :: a = 1.0_RWP
   !> Type of boundary conditions. Currently only PERIODIC is supported.
   INTEGER(KIND=IWP),PARAMETER :: bcs = PERIODIC
   !> Time frequency with which to write data files.
   REAL(KIND=RWP),PARAMETER :: writePeriod = 0.2_RWP
   !> Time frequency with which to print information to the user.
   REAL(KIND=RWP),PARAMETER :: printPeriod = 0.2_RWP
   !> Whether or not to report error analysis. The user might not want to report
   !! errors if the end time is not an integer multiple of the wave speed.
   LOGICAL,PARAMETER :: reportErrors = .TRUE.

   ! Variables required for simulation.
   !
   !> Number of ghost cells at the boundary of the domain.
   INTEGER(KIND=IWP) :: ngc
   !> The computational grid. Note how for a given cell 'i', x(i) is the x
   !! location of the right cell face.
   REAL(KIND=RWP),DIMENSION(0:n) :: x
   !> The locations of the cell centers.
   REAL(KIND=RWP),DIMENSION(1:n) :: xc
   !> The grid spacing.
   REAL(KIND=RWP) :: dx
   !> Solution array at the start and end of the time step.
   REAL(KIND=RWP),DIMENSION(:),ALLOCATABLE :: us
   !> Solution array for the intermediate RK steps.
   REAL(KIND=RWP),DIMENSION(:),ALLOCATABLE :: ui
   !> Array to store the initial condition for error checking.
   REAL(KIND=RWP),DIMENSION(:),ALLOCATABLE :: ic
   !> Array for the RHS of the ODE system.
   REAL(KIND=RWP),DIMENSION(1:n) :: Lu
   !> Time step.
   REAL(KIND=RWP) :: dt
   !> Current simulation time.
   REAL(KIND=RWP) :: t
   !> Left most index in arrays allocated with ghost layers.
   INTEGER(KIND=IWP) :: iL
   !> Right most index in arrays allocated with ghost layers.
   INTEGER(KIND=IWP) :: iR
   !> Current simulation step.
   INTEGER(KIND=IWP) :: nadv = 0_IWP
   !> Whether or not to write data this step.
   LOGICAL :: writeBool = .FALSE.
   !> Whether or not to print information to the user this step.
   LOGICAL :: printBool = .FALSE.
   !> Time after which a data file will be written.
   REAL(KIND=RWP) :: writeTime = writePeriod
   !> Time after which information will be printed to the user.
   REAL(KIND=RWP) :: printTime = printPeriod

   ! Extraneous variables.
   !
   !> Looping index.
   INTEGER(KIND=IWP) :: i
   !> Indices used to set the box wave initial conditions.
   INTEGER(KIND=IWP) :: n1, n2
   !> Variables for order of convergence checking.
   REAL(KIND=RWP) :: l2Err, lInfErr, diff
   !> Logic to determine when to loop and when to stop the loop.
   LOGICAL :: loopBool, exitNextStep

   ! Determine the number of ghost layers based on the spatial order.
   !
   ! NOTE: this will have to change when non-periodic BCs are implemented.
   SELECT CASE (k)
      CASE (2_IWP)
         ngc = 1_IWP
         iL = 0_IWP
         iR = n + ngc
      CASE (3_IWP)
         ngc = 2_IWP
         iL = -1_IWP
         iR = n + ngc
      CASE DEFAULT
         WRITE(*,*) 'Order of reconstruction must be 2 or 3. Halting'
         STOP
   END SELECT

   ! Set the reconstruction order/parameters in the flux module.
   CALL SetReconstructionOrder(n, k)

   ! Allocate solution arrays at various time steps.
   ALLOCATE(us(iL:iR))
   ALLOCATE(ui(iL:iR))
   us(:) = 0.0_RWP
   ui(:) = 0.0_RWP

   ! Fill in the computational grids
   dx = L/REAL(n, RWP)
   DO i = 0, n-1
      x(i) = xstart + REAL(i, RWP)*dx
      xc(i+1) = dx/2.0_RWP + REAL(i, RWP)*dx
   END DO
   x(n) = xstart + L

   ! Form the initial conditions.
   SELECT CASE (ics)
      CASE (SINUSOIDAL)
         ! The initial condition is the integral average of sin(2*pi*x/L).
         DO i = 1, n
            us(i) = L/(2.0_RWP*PI*dx)*(COS(2.0_RWP*PI*x(i-1)/L) - &
                                       COS(2.0_RWP*PI*x(i)/L))
         END DO
      CASE (SQUARE_WAVE)
         ! The initial condition is a square wave.
         n1 = n/3
         n2 = 2*n1
         us(1:n1) = 0.0_RWP
         us(n1+1:n2) = 1.0_RWP
         us(n2+1:n) = 0.0_RWP
      CASE DEFAULT
         WRITE(*,*) 'Invalid IC option. Defaulting to sinusoidal.'
         us(1:n) = SIN(2.0_RWP*PI*xc(1:n)/L)
   END SELECT
   !
   ! Write out the initial conditions.
   CALL WriteData(xc, us, n, ngc, t, nadv)
   !
   ! Store the initial condition for error checking.
   ALLOCATE(ic(1:n))
   ic(:) = us(1:n)

   ! Determine the first time step.
   IF (dtConst < 0.0_RWP) THEN
      SELECT CASE (sys)
         CASE (ADVECTION)
            ! Here we use the advection speed specified as an input.
            dt = cfl*dx/a
         CASE (BURGERS)
            ! For Burgers' the varaible u is actually the wavespeed, so use it.
            dt = cfl*dx/MAXVAL(ABS(us))
         CASE DEFAULT
            WRITE(*,*) 'Invalid system option. Halting'
            STOP
      END SELECT
   ELSE
      dt = dtConst
   END IF

   ! Print some information to the user.
   WRITE(*,100) '-------------------------------------------------------'
   WRITE(*,100) 'Simple1D: A 1D Finite Volume Code for Conservation Laws'
   WRITE(*,100) '-------------------------------------------------------'
   WRITE(*,100) ''
   SELECT CASE (ADVECTION)
      CASE (0)
         WRITE(*,110) 'System:', 'scalar advection'
         WRITE(*,120) 'Wave speed:', a
      CASE (1)
         WRITE(*,110) 'System:', 'Burgers equation'
   END SELECT
   WRITE(*,130) 'Poly. order:', k
   WRITE(*,120) 'Domain length:', L
   WRITE(*,130) 'Num. of cells:', n
   IF (dtConst < 0.0_RWP) THEN
      WRITE(*,120) 'CFL number:', cfl
   ELSE
      WRITE(*,120) 'Fixed dt:', dtConst
   END IF
   IF (nEnd < 0_IWP) THEN
      WRITE(*,120) 'End time:', tend
   ELSE
      WRITE(*,130) 'End iteration:', nEnd
   END IF
   SELECT CASE (bcs)
      CASE (PERIODIC)
         WRITE(*,110) 'BCs:', 'periodic'
   END SELECT
   SELECT CASE (ics)
      CASE (SINUSOIDAL)
         WRITE(*,110) 'ICs:', 'sinusoidal'
      CASE (SQUARE_WAVE)
         WRITE(*,110) 'ICs:', 'Riemann'
   END SELECT
   WRITE(*,100) ''
   100 FORMAT (A)
   110 FORMAT (A,T17,A)
   120 FORMAT (A,T16,ES15.8)
   130 FORMAT (A,T17,I10.10)

   ! Enter the main time loop, in which a 3rd order TVD RK scheme is used.
   t = 0.0_RWP
   nadv = 0_IWP
   loopBool = .TRUE.
   exitNextStep = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Fill in the boundary conditions.
      CALL BoundaryConditions(us, n, ngc, bcs)
      !
      ! Calculate RHS with the data at the start of the time step.
      CALL Flux(us, a, n, k, ngc, dx, Lu)
      !
      ! Determine the first intermediate stage in the RK scheme.
      ui(1:n) = 0.0_RWP
      ui(1:n) = us(1:n) + dt*Lu(1:n)

      ! Fill in the boundary conditions.
      CALL BoundaryConditions(ui, n, ngc, bcs)
      !
      ! Calculate the RHS using the first intermediate stage.
      CALL Flux(ui, a, n, k, ngc, dx, Lu)
      !
      ! Determine the second intermediate stage in the RK scheme.
      ui(1:n) = 0.75_RWP*us(1:n) + 0.25_RWP*ui(1:n) + dt*0.25_RWP*Lu(1:n)

      ! Fill in the boundary conditions.
      CALL BoundaryConditions(ui, n, ngc, bcs)
      !
      ! Calculate the RHS using the second intermediate stage.
      CALL Flux(ui, a, n, k, ngc, dx, Lu)
      !
      ! Determine u at the next time step.
      us(1:n) = us(1:n)/3.0_RWP + &
                2.0_RWP*ui(1:n)/3.0_RWP + &
                dt*2.0_RWP*Lu(1:n)/3.0_RWP

      ! Check whether or not we will write out the data to file.
      IF (writeBool) THEN
         CALL WriteData(xc, us, n, ngc, t, nadv)
         writeTime = writeTime + writePeriod
         writeBool = .FALSE.
      ELSE
         IF (t > writeTime) THEN
            writeBool = .TRUE.
         END IF
      END IF

      ! Check whether or not we will print information to the user.
      IF (printBool .OR. (t == tend))  THEN
         printTime = printTime + printPeriod
         printBool = .FALSE.
         WRITE(*,20) 'Simulation step number: ', nadv, &
                     '; Simulation time: ', t, &
                     '; Max: ', MAXVAL(us(1:n)), &
                     '; Min: ', MINVAL(us(1:n))
         20 FORMAT (A,I10.10,A,ES15.8,A,ES15.8,A,ES15.8)
      ELSE
         IF (t > printTime) THEN
            printBool = .TRUE.
         END IF
      END IF

      ! We have now gone one full step, so increase the step and time counters.
      nadv = nadv + 1_IWP
      t = t + dt

      ! Calculate the time step for the next iteration, if needed.
      IF (dtConst < 0.0_RWP) THEN
         ! The time step can only change for Burgers' equation.
         IF (sys == BURGERS) THEN
            ! For Burgers' the varaible u is the wavespeed.
            dt = cfl*dx/MAXVAL(ABS(us))
         END IF
      END IF

      ! Evaluate exit conditions.
      IF (nEnd > 0_IWP) THEN
         ! If the number of steps is the desired criterion.
         IF (nadv == nEnd) THEN
            EXIT tloop
         END IF
      ELSE
         ! If the end time is the desired criterion.
         IF (exitNextStep) THEN
            EXIT tloop
         END IF
         !
         ! Adjust the time step so the desired end time is exactly reached.
         IF (t + dt > tend) THEN
            dt = tend - t
            exitNextStep = .TRUE.
         END IF
      END IF
   END DO tloop

   ! Write out the final simulation step and time to the user.
   WRITE(*,100) ''
   WRITE(*,130) 'Final step:', nadv
   WRITE(*,120) 'Final time:', t

   ! Write out the final solution.
   CALL WriteData(xc, us, n, ngc, t, nadv)

   ! Calculate the errors.
   IF (reportErrors) THEN
      l2Err = 0.0_RWP
      lInfErr = 0.0_RWP
      SELECT CASE (ics)
         CASE (SINUSOIDAL)
            DO i = 1, n
               diff = us(i) - ic(i)
               l2Err = l2Err + diff**2
               IF (ABS(diff) > lInfErr) THEN
                  lInfErr = ABS(diff)
               END IF
            END DO
         CASE (SQUARE_WAVE)
            DO i = 1, n
               IF ((i >= 1_IWP) .AND. (i <= n1)) THEN
                  l2Err = l2Err + us(i)**2
               ELSE IF ((i > n1) .AND. (i <= n2)) THEN
                  l2Err = l2Err + (us(i) - 1.0_RWP)**2
               ELSE
                  l2Err = l2Err + us(i)**2
               END IF
            END DO
      END SELECT
      l2Err = SQRT(l2Err*dx)
      WRITE(*,100) ''
      WRITE(*,120) 'L2 Error:', l2Err
      WRITE(*,120) 'LInf Error:', lInfErr
   END IF

   ! Deallocate memory.
   DEALLOCATE(us)
   DEALLOCATE(ui)
   DEALLOCATE(ic)
   CALL FluxFinalize()

CONTAINS

   !> Routine to write the solution to file.
   !!
   !> @param[in] xc Locations of cell centers.
   !> @param[in] us Solution to be written.
   !> @param[in] n Number of cells in the domain.
   !> @param[in] ngc Number of ghost layers in the domain.
   !> @param[in] t Current simulation time.
   !> @param[in] nadv Current simulation step.
   SUBROUTINE WriteData(xc, us, n, ngc, t, nadv)
      IMPLICIT NONE
      ! Calling arguments.
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: xc
      REAL(KIND=RWP),DIMENSION(1-ngc:n+ngc),INTENT(IN) :: us
      INTEGER(KIND=IWP),INTENT(IN) :: n, ngc, nadv
      REAL(KIND=RWP),INTENT(IN) :: t
      ! Local variables.
      ! Output file name.
      CHARACTER(LEN=256) :: fname
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Form the output file name.
      WRITE(fname,10) nadv, '.dat'
      10 FORMAT (I12.12,A)

      ! Loop over the cells and record the data to file.
      OPEN(UNIT=100,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      DO i = 1, n
         WRITE(100,20) xc(i), us(i)
         20 FORMAT (ES15.8,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteData

END PROGRAM Simple1D

