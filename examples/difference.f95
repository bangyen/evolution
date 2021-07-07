! Integro diff eq. 13 oct 2014 (Reformatted)
! http://www.uprh.edu/rbaretti/IntegroDE13oct2014.htm

! declarations: types, data, and print format
      implicit real*8 (a-h, o-z)
      data nstep, L /2000, 4.d0/
  100 format('x, phi, phiex, ratio = ' 4(3x, e10.3))

! initial values:
      phiex(x) = (1.d0 / 2.d0)    &
                 * dexp(-x)       &
                 * dsin(2.d0 * x)
      dx       = L / dfloat(nstep)
      kp       = int(dfloat(nstep) / 40.d0)
      kount    = kp
      phi0     = 0.d0
      uprime0  = 1.d0
      phi1     = phi0 + dx * uprime0

! print initial values
      print 100,               &
            0.d0,              &
            phi0,              &
            phiex(0.d0),       &
            phi0 / phiex(0.d0)

! for loop
      do i = 2, nstep
          x      = dx * dfloat(i)
          dphidx = (phi1 - phi0) / dx

! because i >= 2, this block is never entered,
! making the if statement is unnecessary
          if (i .eq. 1) then
              phi2 = 2.d0 * phi1         &
                     - phi0              &
                     + dx ** 2           &
                     * (- 2.d0 * dphidx  &
                        - 5.d0 * phi1    &
                        + 1.d0 / dx)
          else
              phi2 = 2.d0 * phi1         &
                     - phi0              &
                     + dx ** 2           &
                     * (- 2.d0 * dphidx  &
                        - 5.d0 * phi1)
          end if
! every kp number of iterations, print
          if (i .eq. kount) then
              print 100,            &
                    x,              &
                    phi2,           &
                    phiex(x),       &
                    phi2 / phiex(x)
              kount = kount + kp
          end if

          phi0 = phi1
          phi1 = phi2
      end do
      stop
      end


! approximation for dirac delta function
! 1/dx near zero, 0 elsewhere
      function delta(x, dx)
          implicit real*8 (a-h, o-z)
          if (x .le. dx) then
              delta = 1.d0 / dx
          else
              delta = 0.d0
          end if
          return
      end
