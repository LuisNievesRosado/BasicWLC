subroutine betacalc(XIR,N,Beta)
    double precision betamat(504,2) ! Matrix of beta values
    double precision Beta           ! Hydrodynamic interaction factor
    double precision XIR            ! Friction of a bead
    double precision L              ! Length of a polymer

    integer N                       ! Number of beads in polymer
    integer i,j                     ! Loop counters

    L = XIR*N

    open(unit = 11, file="input/beta")
    do i=1,504
       read(11,*) betamat(i,1:2)
    end do
    close(11)

    if (2*betamat(504,1).LT.L) then
       Beta = betamat(504,2)
    else if (2*betamat(1,1).GE.L) then
       print*, "Values of beta are not available for chains shorter than 4 lp"
    else
       do j=1,504
          if (2*betamat(j,1).EQ.L) then
             Beta = betamat(j,2)
          else if (j.NE.1.AND.j.NE.504.AND.2*betamat(j-1,1).LT.L.AND.2*betamat(j+1,1).GT.L) then
             Beta = (betamat(j-1,2)+betamat(j+1,2))/2d0
          end if
       end do
    end if

end subroutine betacalc
