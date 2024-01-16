!=====================================================================================
! tonhonr@usp.br
!
! This module contains the most fundamental observables of pure gauge theory
! We implement the planar Wilson loops and the Polyacov loop
!
!=====================================================================================
module measures
    use lattice
    use mathtools
    implicit none
    real(kind=r8), allocatable, dimension(:,:) :: w, wp!, wr
    real(kind=r8), allocatable, dimension(:) :: vqq, vqqp, chi, chip!, chip
    contains
    subroutine smearing(u,alpha,nsmearing)
        complex(kind=r8), dimension(nr,nr,nr,nt,4,nc,nc), intent(inout) :: u
        complex(kind=r8) :: uprime(nr,nr,nr,nt,4,nc,nc), umini(nc,nc)
        real(kind=r8), intent(in) :: alpha
        integer(kind=i4) :: mi, x, y, z, t, nsmearing, i
        ! number of smearing sweeps
        do i=1,nsmearing
            uprime=cz
            ! access each lattice link
            do x=1,nr
                do y=1,nr
                    do z=1,nr
                        do t=1,nt
                            do mi=1,4
                                ! do the smearing sweep
                                umini=stample(u,x,y,z,t,mi)
                                uprime(x,y,z,t,mi,:,:)=(1.0_r8-alpha)*u(x,y,z,t,mi,:,:)+alpha*dagger(umini)/6.0_r8
                                call keapinsuN(uprime(x,y,z,t,mi,:,:),nc)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
            u=uprime
        enddo
    endsubroutine smearing

    subroutine measurewilson(u,a,b,w)
        complex(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,nc,nc)
        real(kind=r8), intent(out) :: w(a,b)
        real(kind=r8) :: sum
        integer(kind=i4) :: a,b, i, j, n(4), e1, e2, e3, e4, mu, nu
        ! we compute the loops w(a,b) and save the result
        do i=1,a
            do j=1,b
                sum=0.0_r8
                ! compute the mean wilson loop
                do e1=1,nr
                do e2=1,nr
                do e3=1,nr
                do e4=1,nt
                    n=(/e1,e2,e3,e4/)
                    do mu=1,3
                    do nu=mu+1,4
                        sum=sum+wilsonplanar(u,n,mu,nu,i,j)
                        ! the lattice is simmetric, so we
                        ! must consider the case a->b and b->a
                        ! now we interchange a and b
                        if(i.ne.j)then 
                            sum=sum+wilsonplanar(u,n,mu,nu,j,i)                    
                        endif
                        ! the loop for this site is calculated!
                    enddo
                    enddo
                enddo
                enddo
                enddo
                enddo
            
                if(i.eq.j)then
                    w(i,j)=sum/(6.0_r8*nt*nr**3)    ! quadratic loops
                else
                    w(i,j)=sum/(12_r8*nt*nr**3) ! rectangular loops
                endif
            enddo
        enddo
    endsubroutine measurewilson

    ! function that calculate te wilson loop of size a,b
    function wilsonplanar(u,n,mu,nu,a,b)
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc)
        complex(kind=r8), dimension(nc,nc) :: umunu
        real(kind=r8) :: wilsonplanar
        integer(kind=i4), dimension(4) :: n
        integer(kind=i4) :: mu, nu, a, b, i
        
            umunu=ident(nc)
            !n=(/x,y,z,t/)
            ! computing the lines
            do i=1,a
            umunu=matmul(umunu,u(n(1),n(2),n(3),n(4),mu,:,:))
            n(mu)=n(mu)+1
            if(n(mu).gt.nr) n(mu)=1
            enddo
            do i=1,b
                umunu=matmul(umunu,u(n(1),n(2),n(3),n(4),nu,:,:))
                n(nu)=n(nu)+1
                if(n(nu).gt.nr) n(nu)=1
            enddo
            do i=1,a
                n(mu)=n(mu)-1
                if(n(mu).lt.1) n(mu)=nr
                umunu=matmul(umunu,conjg(transpose(u(n(1),n(2),n(3),n(4),mu,:,:))))
            enddo
            do i=1,b
                n(nu)=n(nu)-1
                if(n(nu).lt.1) n(nu)=nr
                umunu=matmul(umunu,conjg(transpose(u(n(1),n(2),n(3),n(4),nu,:,:))))
            enddo
            wilsonplanar=dreal(ctrace(umunu,nc))/nc
    endfunction wilsonplanar

    ! subroutine that compute the creutz ratios
    subroutine creutz(chi,a,b,w)
        real(kind=r8) :: w(a,b), chi(a)
        integer(kind=i4) :: a, b, i
        chi=0.0_r8
        do i=1,a
            if(i.eq.1)then
                chi(1)=-log(w(1,1))
            else
                chi(i)=-log((w(i,i)*w(i-1,i-1))/(w(i-1,i)*w(i,i-1)))
            endif
        enddo
    endsubroutine creutz

    ! subroutine tat computes the qq potential on the lattice
    subroutine qqpot(v,w,a,b)
        real(kind=r8) :: v(a), w(a,b)
        integer(kind=i4) :: a, b, i, j
        v=0.0_r8
        do i=1,a
        do j=1,b-1
            v(i)=v(i)-log(w(i,j+1)/w(i,j))
        enddo
        enddo
        v=v/(b-1)
    endsubroutine qqpot
endmodule measures