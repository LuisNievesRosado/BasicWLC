! This subroutine computes the Rotne-Prager-Yamakawa diffusivity
! tensor for beads moving under hydrodynamic interactions.
! Author: Luis Nieves
! Last Edit: August 9, 2016

subroutine diffcalc(R,NT,XIR,a,Beta,D)

implicit none
double precision D(NT,NT,3,3) ! RPY Diffusion tensor
double precision R(NT,3)      ! Bead position
double precision XIR          ! Bead friction factor
double precision Beta         ! Hydrodynamic Interaction factor

double precision a            ! Hydrodynamic radius of bead
double precision rx,ry,rz,rij ! Relative position components

integer NT                    ! Number of beads
integer I,J                   ! Loop counters


do I = 1,NT

    D(I,I,1,1) = Beta/XIR
	D(I,I,2,2) = Beta/XIR
	D(I,I,3,3) = Beta/XIR
	
    do J = 1,I-1
	
	    rx(1) = R(I,1)-R(J,1)
		ry(2) = R(I,2)-R(J,2)
		rx(3) = R(I,3)-R(J,3)
		rij = sqrt(rx^2+ry^2+rx^2)
		
		if (rij < 2*a) then
		
	        D(I,J,1,1) = rx^2/8/a^2/rij + 8/6/a - 3*rij/8/a^2
	        D(I,J,2,1) = rx*ry/8/a^2/rij
		    D(I,J,2,2) = ry^2/8/a^2/rij + 8/6/a - 3*rij/8/a^2
			D(I,J,3,1) = rx*rz/8/a^2/rij
		    D(I,J,3,2) = ry*rz/8/a^2/rij
		    D(I,J,3,3) = rz^2/8/a^2/rij + 8/6/a - 3*rij/8/a^2
			             
	    else             
                         
			D(I,J,1,1) = rx^2/rij^3*(1-2*a^2/rij^2) + 1/rij + 2*a^2/3/rij^3
	        D(I,J,2,1) = rx*ry/rij^3*(1-2*a^2/rij^2)
		    D(I,J,2,2) = ry^2/rij^3*(1-2*a^2/rij^2) + 1/rij + 2*a^2/3/rij^3
			D(I,J,3,1) = rx*rz/rij^3*(1-2*a^2/rij^2)
		    D(I,J,3,2) = ry*rz/rij^3*(1-2*a^2/rij^2)
		    D(I,J,3,3) = rz^2/rij^3*(1-2*a^2/rij^2) + 1/rij + 2*a^2/3/rij^3
			
    	end if 
	end do
end do	

end subroutine diffcalc
