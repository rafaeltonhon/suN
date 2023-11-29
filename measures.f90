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
    real(kind=r8), allocatable :: w(:,:)
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

    subroutine measurewilson(u,r,t,w)
        complex(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,nc,nc)
        real(kind=r8), intent(inout) :: w(r,t)
        integer(kind=i4) :: r, t, i, j
        ! we compute the loops w(a,b) and save the result
        w=0.0_r8
        do i=1,r
            do j=1,t
                w(i,j)=wilsonplanar(u,i,j)
            enddo
        enddo
    endsubroutine measurewilson

    ! function that calculate te wilson loop of size a,b
    function wilsonplanar(u,a,b)
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc)
        complex(kind=r8), dimension(nc,nc) :: unmu, unpmunu, unpnumu, unnu, umunu
        real(kind=r8) :: sum, wilsonplanar
        integer(kind=i4), dimension(4) :: n, m, npmu, npnu
        integer(kind=i4) :: x, y, z, t, mu, nu, a, b, i
        sum=0.0_r8
        do x=1,nr
            do y=1,nr
                do z=1,nr
                    do t=1,nt
                        do mu=1,3
                            do nu=mu+1,4
                                umunu=ident(nc)
                                n=(/x,y,z,t/)
                                m=(/x,y,z,t/)
                                unmu=ident(nc)
                                unpmunu=ident(nc)
                                unpnumu=ident(nc)
                                unnu=ident(nc)


                                ! computing the lines
                                do i=1,a
                                    !unmu=matmul(unmu,u(n(1),n(2),n(3),n(4),mu,:,:))
                                    umunu=matmul(umunu,u(n(1),n(2),n(3),n(4),mu,:,:))
                                    n(mu)=n(mu)+1
                                    if(n(mu).gt.nr) n(mu)=1
                                enddo

                                do i=1,b
                                    !unpmunu=matmul(unpmunu,u(n(1),n(2),n(3),n(4),nu,:,:))
                                    umunu=matmul(umunu,u(n(1),n(2),n(3),n(4),nu,:,:))
                                    n(nu)=n(nu)+1
                                    if(n(nu).gt.nr) n(nu)=1
                                enddo

                                do i=1,a
                                    n(mu)=n(mu)-1
                                    if(n(mu).lt.1) n(mu)=nr
                                !   unpnumu=matmul(unpnumu,u(n(1),n(2),n(3),n(4),mu,:,:))
                                umunu=matmul(umunu,conjg(transpose(u(n(1),n(2),n(3),n(4),mu,:,:))))
                                enddo
                                !umunu=matmul(umunu,conjg(transpose(unpnumu)))
                                
                                do i=1,b
                                    n(nu)=n(nu)-1
                                    if(n(nu).lt.1) n(nu)=nr
                                    umunu=matmul(umunu,conjg(transpose(u(n(1),n(2),n(3),n(4),nu,:,:))))
                                    !unnu=matmul(unnu,u(n(1),n(2),n(3),n(4),nu,:,:))
                                enddo
                                sum=sum+dreal(ctrace(umunu,nc))/(nc)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        if(a.eq.b)then
            wilsonplanar=sum/(6.0_r8*nt*nr**3)    ! quadratic loops
        else
            wilsonplanar=sum/(12_r8*nt*nr**3) ! rectangular loops
        endif
    endfunction wilsonplanar

    ! function that calculates the polyakov lop at a position m
    function polyakovl(x,y,z,nt,u)
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc), pl(nc,nc)
        complex(kind=r8) :: polyakovl
        integer(kind=i4) :: x, y, z, t, nt
        pl=ident(nc)
        do t=1,nt
            pl=matmul(pl,u(x,y,z,t,4,:,:))
        enddo
        polyakovl=(ctrace(pl,nc))!pl(1,1)+pl(2,2))
    endfunction polyakovl

    ! function that retunrns a matrix with the polyakov loop on each side of the space
    function latticepolyakov(nr,nt,u)
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc)
        complex(kind=r8) :: plat(nr,nr,nr), latticepolyakov(nr,nr,nr)
        integer(kind=i4) :: x, y, z, nr, nt
        do x=1,nr
            do y=1,nr
                do z=1,nr
                    plat(x,y,z)=polyakovl(x,y,z,nt,u)
                enddo
            enddo
        enddo
        latticepolyakov=plat
    endfunction latticepolyakov

    subroutine polyakovcorrelator(pmpn,nr)
        real(kind=r8), intent(inout) :: pmpn(nr)
        complex(kind=r8) :: polyakov(nr,nr,nr)
        integer(kind=i4) :: nr, x, y, z, r, mi, m(3), n(3)
        ! computes the polyakov loops on each site of the lattice
        polyakov=latticepolyakov(nr,nt,u)
        pmpn=0.0_r8
        ! no we do the mean value of the correlator
        do x=1,nr
            do y=1,nr
                do z=1,nr
                    m=(/x,y,z/)
                    ! we vary the distance
                    do r=1,nr
                        ! we varry the direction
                        do mi=1,3
                            n=m
                            n(mi)=n(mi)+r
                            if(n(mi).gt.nr) n(mi)=n(mi)-nr
                            pmpn(r)=pmpn(r)+dreal(polyakov(x,y,z)*conjg(polyakov(n(1),n(2),n(3))))
                        enddo
                    enddo
                enddo
            enddo
        enddo
        pmpn=pmpn/(2.0_r8*nr**3.0_r8)
    endsubroutine polyakovcorrelator
endmodule measures