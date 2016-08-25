! This subroutine computes the brownian mobility tensor
! using a Cholesky decomposition.

subroutine brownchol(D,NT,sigB)

implicit none
double precision D(NT,NT,3,3)    ! RPY Mobility tensor
double precision sigB(NT,NT,3,3) ! Brownian Mobility tensor
double precision sigsum          ! Summation variable

integer NT                       ! Number of beads
integer I,Io,Ii                  ! Loop counters       
integer J,Jo,Ji                  ! Loop counters
integer K,Ko,Ki                  ! Loop counters

do I = 1,3*NT

     Io = I/3 + 1
     Ii = 3 - 3*Io + I 
 
     do J = I,3*NT

        Jo = J/3 + 1
        Ji = 3 - 3*Jo + J
 
        if (Ii > Ji) then
           sigsum = D(Io,Jo,Ii,Ji)
        else
           sigsum = D(Io,Jo,Ji,Ii) 
        end if
 
        do K = I-1,1,-1
           Ko = K/3 + 1
           Ki = 3 - 3*Ko + K 
           sigsum = sigsum - sigB(Io,Ko,Ii,Ki)*sigB(Jo,Ko,Ji,Ki)
        end do
 
        if (I.EQ.J) then
           sigB(Io,Io,Ii,Ii) = sqrt(sigsum)
        else 
           sigB(Jo,Io,Ji,Ii) = sigsum/sigB(Io,Io,Ii,Ii)

        end if
     end do 
end do

end subroutine brownchol
