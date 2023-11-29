module lattice
    implicit none

    ! defining our types
    integer, parameter :: r8=selected_real_kind(15,8)
    integer, parameter :: i4=selected_int_kind(8)

    ! defining our parameters
    complex(kind=r8), parameter :: cr=(1.0_r8,0.0_r8)
    complex(kind=r8), parameter :: ci=(0.0_r8,1.0_r8)
    complex(kind=r8), parameter :: cz=(0.0_r8,0.0_r8)
    real(kind=r8), parameter :: pi=dacos(-1.0_r8)
    !integer(kind=i4), parameter :: nt=16
    !integer(kind=i4), parameter :: nr=16
    !integer(kind=i4), parameter :: nc=2

    ! defining the needed matrices
    complex(kind=r8), allocatable, dimension(:,:,:,:,:,:,:) :: u

    ! defining our variables
    real(kind=r8) :: x1, x2, x3, x4, beta, epsilon
    integer(kind=i4) :: e1, e2, e3, e4, mi, ni
    integer(kind=i4) :: nconf, nr, nt, nc

    public :: nr, nt, x1, x2, x3, x4, beta
    public :: e1, e2, e3, e4, mi, ni
    contains
     ! function that calculates the sum of the stamples
    function stample(u,e1,e2,e3,e4,mi)
        complex(kind=r8), dimension(nc,nc) :: stample, uni_pmi, umi_pnid, uni_nd, p1, p2
        complex(kind=r8), dimension(nc,nc) :: uni_pmimnid, umi_mnid, uni_mni, temp
        integer(kind=i4), dimension(1:4) :: n, npmi, nmmi, npni, nmni, npmimni
        complex(kind=r8) :: u(nr,nr,nr,nt,4,nc,nc)
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni

        temp=0.0_r8
        n=(/e1,e2,e3,e4/)
        do ni=1,4
            if(ni.ne.mi)then
                call neighbors(n,npmi,nmmi,npni,nmni,npmimni,mi,ni)

                uni_pmi=u(npmi(1),npmi(2),npmi(3),npmi(4),ni,:,:)
                umi_pnid=conjg(transpose(u(npni(1),npni(2),npni(3),npni(4),mi,:,:)))
                uni_nd=conjg(transpose(u(n(1),n(2),n(3),n(4),ni,:,:)))

                uni_pmimnid=conjg(transpose(u(npmimni(1),npmimni(2),npmimni(3),npmimni(4),ni,:,:)))
                umi_mnid=conjg(transpose(u(nmni(1),nmni(2),nmni(3),nmni(4),mi,:,:)))
                uni_mni=u(nmni(1),nmni(2),nmni(3),nmni(4),ni,:,:)

                temp=temp+matmul(uni_pmi,matmul(umi_pnid,uni_nd))
                temp=temp+matmul(uni_pmimnid,matmul(umi_mnid,uni_mni))
            endif
        enddo

        stample=temp
    endfunction stample

        ! subroutine that return us the six neighbors
    subroutine neighbors(n,npmi,nmmi,npni,nmni,npmimni,mi,ni)
        integer(kind=i4), dimension(1:4) :: n, npmi, nmmi, npni, nmni, npmimni
        integer(kind=i4) :: mi, ni
        npmi=n
        nmmi=n
        npni=n
        nmni=n

        ! calculating the vectors with the neighbors positionr
        npmi(mi)=npmi(mi)+1
        nmmi(mi)=nmmi(mi)-1
        npni(ni)=npni(ni)+1
        nmni(ni)=nmni(ni)-1

        npmimni=npmi
        npmimni(ni)=npmimni(ni)-1

        ! applyng the periodic bounduary conditionr
        if((npmi(mi).gt.nr).or.(npmi(mi).gt.nt))then
            npmi(mi)=1
            npmimni(mi)=1
        endif
        if((npni(ni).gt.nr).or.(npni(ni).gt.nt))then
            npni(ni)=1
        endif
        if(nmmi(mi).lt.1)then
            if(mi.eq.4)then
                nmmi(mi)=nt
            else
                nmmi(mi)=nr
            endif
        endif
        if(nmni(ni).lt.1)then
            if(ni.eq.4)then
                nmni(ni)=nt
                npmimni(ni)=nt
            else
                nmni(ni)=nr
                npmimni(ni)=nr
            endif
        endif
    endsubroutine neighbors
endmodule lattice