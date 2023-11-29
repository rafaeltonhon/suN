module mathtools
    use lattice
    implicit none
    contains
    ! returns an nc x nc identity matrix
    function ident(nc)
        complex(kind=r8), dimension(nc,nc) :: ident, temp
        integer(kind=i4) :: i, nc
        temp=cz
        do i=1,nc
            temp(i,i)=cr
        enddo
        ident=temp
    endfunction ident
    
    ! function that returns the trace of a complex matrix
    function ctrace(m,nc)
        complex(kind=r8) :: m(nc,nc), t, ctrace
        integer(kind=i4) :: i, nc
        t=cz
        do i=1,nc
            t=t+m(i,i)
        enddo
        ctrace=t
    endfunction ctrace

    ! function that returns the hermitian of a matrix
    function dagger(m)
        complex(kind=r8), dimension(nc,nc) :: dagger, m
        dagger=conjg(transpose(m))
    endfunction dagger

    ! function tht returns the inner product
    function inner(u,v)
        complex(kind=r8) :: inner, temp, u(3), v(3)
        temp=u(1)*conjg(v(1))+u(2)*conjg(v(2))+u(3)*conjg(v(3))
        inner=temp
    endfunction inner

    ! funtion that returns the cross product of two vectors
    function crossprod(a,b)
        complex(kind=r8), dimension(3) :: crossprod, a, b, c
        integer(kind=i4) :: i, j, k
        c=cz
        do i=1,3
            do j=1,3
                do k=1,3
                    c(k)=c(k)+a(i)*b(j)*eps(i,j,k)
                enddo
            enddo
        enddo
        crossprod=c
    endfunction crossprod

    !
    function eps(i,j,k)
        integer(kind=i4) :: eps, i, j, k
        if((i.eq.j).or.(i.eq.k).or.(j.eq.k))then
            eps=0
        elseif((i.eq.1).and.(j.eq.2).and.(k.eq.3))then
            eps=1
        elseif((i.eq.3).and.(j.eq.1).and.(k.eq.2))then
            eps=1
        elseif((i.eq.2).and.(j.eq.3).and.(k.eq.1))then
            eps=1
        else
            eps=-1
        endif
    endfunction eps

    ! 2d Gram-Schmidt process
    subroutine gramschmidt2d(a,b)
        complex(kind=r8), dimension(3), intent(inout) :: a, b
        real(kind=r8) :: moda, modb
        ! normalize a
        moda=((inner(a,a)))!**0.5_r8
        a=a/dsqrt(moda)
        ! make b orthonormal to a
        b=b-dreal(inner(b,a))*a
        modb=(dreal(inner(b,b)))!**0.5_r8
        b=b/dsqrt(modb)
    endsubroutine gramschmidt2d

    ! subroutine that normilize the link variable
    subroutine keapinsuN(ul,nc)
        complex(kind=r8), intent(inout) :: ul(nc,nc)
        complex(kind=r8) :: u(3), v(3), ucv(3)
        real(kind=r8) :: detul
        integer(kind=i4) :: nc, i
        select case(nc)
        case(2) ! su(2) group
            detul=0.0_r8
            detul=detul+dreal(ul(1,1))**2.0_r8+aimag(ul(1,1))**2.0_r8
            detul=detul+dreal(ul(2,1))**2.0_r8+aimag(ul(2,1))**2.0_r8
            ul=ul/(detul**0.5_r8)
        case(3) ! su(3) group
            u=(/ul(1,1),ul(1,2),ul(1,3)/)
            v=(/ul(2,1),ul(2,2),ul(2,3)/)
            call gramschmidt2d(u,v)
            ucv=crossprod(conjg(u),conjg(v))
            do i=1,3
                ul(1,i)=u(i)
                ul(2,i)=v(i)
                ul(3,i)=ucv(i)
            enddo
        case default
            stop
        end select
    endsubroutine keapinsuN
endmodule mathtools