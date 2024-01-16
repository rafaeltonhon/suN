!=====================================================================================
! tonhonr@usp.br
!
! This module contains routines for center projection and identification of center vortex
! THe gauge fix routines, need for center projection, are on a separed module
!=====================================================================================
module center
    use lattice
    use mathtools
    use measures
    implicit none
    contains
    ! this routine recives a nonprojected configuration u (with gauge alread fixed)
    ! and return us a center projected one, z=sign(a4)I
    ! for now, this only work on su(2)
    subroutine centerprojection(u,z)
        complex(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,nc,nc)
        integer(kind=i4), intent(out) :: z(nr,nr,nr,nt,4)
        integer(kind=i4) :: e1, e2, e3, e4, mi
        ! acess each lattic site
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                ! project the link on the center of the group
                if(real(u(e1,e2,e3,e4,mi,1,1)).gt.0.0_r8) then
                    z(e1,e2,e3,e4,mi)=1
                else
                    z(e1,e2,e3,e4,mi)=-1
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine centerprojection

    ! function that says if an giben paquette has a vortex pircing it
    function paquettepierce(z,n,mi,ni)
        integer(kind=i4) :: z(nr,nr,nr,nt,4)
        integer(kind=i4) :: mi, ni, paquettepierce
        integer(kind=i4), dimension(4) :: n, npmi, npni, v
        call neighbors(n,npmi,v,npni,v,v,mi,ni)
        paquettepierce=z(n(1),n(2),n(3),n(4),mi)*z(npmi(1),npmi(2),npmi(3),npmi(4),ni)*&
                    &z(npni(1),npni(2),npni(3),npni(4),mi)*z(n(1),n(2),n(3),n(4),ni)
    endfunction paquettepierce

    ! function that, given the size of an planar wilson loop at site n
    ! in directions mi, ni tell us how may P-vortices there is inside this loop
    function nvortices(z,n,mi,ni,r,t)
        integer(kind=i4) :: z(nr,nr,nr,nt,4), n(4), m(4), i, j
        integer(kind=i4) :: mi, ni, r, t, nv, nvortices
        nv=0
        m=n
        ! we need to analyse plaquette per plaquette of the loop
        do i=1,r
            do j=1,t
                ! we analyse the plaquettes
                if(paquettepierce(z,n,mi,ni).lt.0) nv=nv+1
                !print*,'plaqu',paquettepierce(z,n,mi,ni)
                n(ni)=n(ni)+1
                if(n(ni).gt.nt) n(ni)=1
            enddo
            n(mi)=n(mi)+1
            if(n(mi).gt.nr) n(mi)=1
            n(ni)=m(ni)
        enddo
        n=m ! do not alter n
        nvortices=nv
    endfunction nvortices

    ! subroutine computes the wilson loops of size r,t
    ! in projected and un projected configurations
    subroutine projectedwilson(u,a,b,wp)
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc)
        real(kind=r8) :: wp(a,b), sum
        integer(kind=i4) :: z(nr,nr,nr,nt,4), n(4), i, j, a, b
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni, lambda
        ! first, we center project the links
        call centerprojection(u,z)
        !wp=0.0_r8
        do i=1,a
        do j=1,b
            lambda=0
            ! we compute the loops w(a,b)
            sum=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mi=1,3
                do ni=mi+1,4
                    ! this loop has a vortex?
                    if(nvortices(z,n,mi,ni,i,j).gt.0)then
                        sum=sum+wilsonplanar(u,n,mi,ni,i,j)
                        lambda=lambda+1
                    endif
                    ! there is no distinction between the directions
                    ! then, we must do i->j and j->i if i!=j
                    if(i.ne.j)then
                        ! this loop has a vortex?
                        n=(/e1,e2,e3,e4/)
                        if(nvortices(z,n,mi,ni,j,i).gt.0)then
                            sum=sum+wilsonplanar(u,n,mi,ni,j,i)
                            lambda=lambda+1
                        endif
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            ! we have computed the loops, letss get away with the result
            if (i.eq.j)then
                wp(i,j)=sum/lambda!(6.0_r8*nr**3.0_r8*nt)
            else
                wp(i,j)=sum/lambda!(12.0_r8*nr**3.0_r8*nt)
            endif
        enddo
        enddo
    endsubroutine projectedwilson
endmodule center