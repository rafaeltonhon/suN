module thermalize
    use lattice
    use mathtools
    implicit none
    contains
    ! function that creates a random matrix x
    function xmatrix(epsilon)
        complex(kind=r8), dimension(2,2) :: xmatrix, temp
        real(kind=r8) :: r(1:3), r4, rn, epsilon
        integer(kind=i4) :: i
        rn=0.0_r8
        call random_number(r4)
        r4=r4-0.5_r8
        do i=1,3
            call random_number(r(i))
            r(i)=r(i)-0.5_r8
            rn=rn+r(i)**2.0_r8
        enddo
        r=r*epsilon/dsqrt(rn)
        if(r4.gt.0.0_r8) r4=(1.0_r8-epsilon**2.0_r8)**0.5_r8
        if(r4.lt.0.0_r8) r4=-(1.0_r8-epsilon**2.0_r8)**0.5_r8

        temp(1,1)=r4*cr+r(3)*ci
        temp(1,2)=r(2)*cr+r(1)*ci
        temp(2,1)=-conjg(temp(1,2))
        temp(2,2)=conjg(temp(1,1))
        xmatrix=temp
    endfunction xmatrix

    function sunxmatrix(epsilon,nc)
        complex(kind=r8), dimension(nc,nc) :: sunxmatrix, r, s, t
        complex(kind=r8) :: x(2,2)
        real(kind=r8) :: epsilon
        integer(kind=i4) :: nc

        select case(nc)
        case(2) ! su2 group
            sunxmatrix=xmatrix(epsilon)
        case(3) ! su3 group
            r=cz
            s=cz
            t=cz

            ! for r
            x=xmatrix(epsilon)
            r(1,1)=x(1,1)
            r(1,2)=x(1,2)
            r(2,1)=x(2,1)
            r(2,2)=x(2,2)
            r(3,3)=cr

            ! for s
            x=xmatrix(epsilon)
            s(1,1)=x(1,1)
            s(1,3)=x(1,2)
            s(3,1)=x(2,1)
            s(3,3)=x(2,2)
            s(2,2)=cr

            ! for t
            x=xmatrix(epsilon)
            t(2,2)=x(1,1)
            t(2,3)=x(1,2)
            t(3,2)=x(2,1)
            t(3,3)=x(2,2)
            t(1,1)=cr

            r=matmul(r,s)
            r=matmul(r,t)
            sunxmatrix=r
        case default
            stop
        end select
    endfunction

    ! subroutne that initialize the gluon matrix and calculates the initial action
    subroutine initialize(u,action,init)
        complex(kind=r8), intent(inout):: u(nr,nr,nr,nt,4,nc,nc)
        complex(kind=r8),dimension(nc,nc) :: plaquette, uni_pmi, umi_pnid, uni_nd, temp
        integer(kind=i4), dimension(1:4) :: n, npmi, nmmi, npni, nmni, npmimni
        real(kind=r8) :: action, sum, r(4), rn, p1
        integer(kind=i4) :: e1, e2, e3, e4, mi, init, i, ni

        sum=0.0_r8
        select case(init)
        case(1) ! we use the cold start
            do e1=1,nr
                do e2=1,nr
                    do e3=1,nr
                        do e4=1,nt
                            do mi=1,4
                                u(e1,e2,e3,e4,mi,:,:)=ident(nc)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        case(2) ! we use the hot start
            ! we initialize this as a random shit
            do e1=1,nr
                do e2=1,nr
                    do e3=1,nr
                        do e4=1,nt
                            do mi=1,4
                                u(e1,e2,e3,e4,mi,:,:)=sunxmatrix(epsilon,nc)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        case default
            stop "novalid initialize option was selected, idiot"
        endselect

        ! now we calculates the action
        do e1=1,nr
            do e2=1,nr
                do e3=1,nr
                    do e4=1,nt
                        do mi=1,4
                            plaquette=cz
                            temp=0.0_r8
                            n=(/e1,e2,e3,e4/)
                            p1=0.0_r8
                            do ni=mi+1,4
                                    p1=p1+1.0_r8
                                    call neighbors(n,npmi,nmmi,npni,nmni,npmimni,mi,ni)
                    
                                    uni_pmi=u(npmi(1),npmi(2),npmi(3),npmi(4),ni,:,:)
                                    umi_pnid=conjg(transpose(u(npni(1),npni(2),npni(3),npni(4),mi,:,:)))
                                    uni_nd=conjg(transpose(u(n(1),n(2),n(3),n(4),ni,:,:)))
                                    temp=temp+matmul(uni_pmi,matmul(umi_pnid,uni_nd))
                            enddo
                            plaquette=matmul(u(e1,e2,e3,e4,mi,:,:),temp)
                            sum=sum+p1*nc*1.0_r8-dreal(ctrace(plaquette,nc))
                        enddo
                    enddo
                enddo
            enddo
        enddo
        action=beta*sum/(nc)
        action=action/(6*nt*nr**3.0_r8)
    endsubroutine initialize

    ! subroutine that does su(2) heat bath hit
    subroutine hbxmatrix(x,alpha)
        complex(kind=r8), dimension(2,2), intent(inout) :: x
        real(kind=r8) :: lambda, alpha, r, r1, r2, r3
        real(kind=r8) :: xi(4), b(3), mb
        integer(kind=i4) :: i, j
        x=cz
        do i=1,500
            call random_number(r)
            call random_number(r1)
            call random_number(r2)
            call random_number(r3)
            r1=1.0_r8-r1
            r2=1.0_r8-r2
            r3=1.0_r8-r3
            lambda=-(dlog(r1)+(dcos(2.0_r8*pi*r2))**2.0_r8*dlog(r3))
            !if(nc==2) lambda=1.0_r8*lambda/(2.0_r8*alpha*beta)
            !if(nc==3) 
            lambda=nc*lambda/(4.0_r8*alpha*beta)
            if(r**2.0_r8.le.(1.0_r8-lambda)) exit
        enddo
        if(i==500) stop 'heatbath xmatrix convergence not achived. line 193'
        xi(4)=1.0_r8-2.0_r8*lambda
        !print*,lambda
        do i=1,500
            mb=0.0_r8
            do j=1,3
                call random_number(b(j))
                b(j)=1.0_r8-2.0_r8*b(j)
                mb=mb+b(j)**2.0_r8
            enddo
                if(mb.le.1.0_r8) exit
        enddo
        if(i==500) stop 'heatbath xmatrix convergence not achived. line 205'

        do i=1,3
            xi(i)=dsqrt(1.0_r8-xi(4)**2.0_r8)*b(i)/dsqrt(mb)
        enddo
        x(1,1)=xi(4)*cr+xi(3)*ci
        x(1,2)=xi(2)*cr+xi(1)*ci
        x(2,2)=xi(4)*cr-xi(3)*ci
        x(2,1)=-xi(2)*cr+xi(1)*ci 
    endsubroutine hbxmatrix
    
    ! subroutine that does one one heat bath hit in su(N) with n>2
    subroutine sunhb(u,a,nc)
        !complex(kind=r8), dimension(nr,nr,nr,nt,4,nc,nc), intent(inout) :: u
        complex(kind=r8), dimension(nc,nc) :: a, r, w, u
        complex(kind=r8) :: w1, w2, v(2,2), y(2,2), x(2,2)
        real(kind=r8) :: detw, alpha
        integer(kind=i4) :: nc
        select case(nc)
        case(2) ! SU(2) group
            alpha=a(1,1)*a(2,2)-a(1,2)*a(2,1)
            alpha=dsqrt(alpha)
            if(dabs(alpha).lt.1e-5)then
                a=xmatrix(0.0_r8)
                alpha=1.0_r8
            endif
            call hbxmatrix(x,alpha)
            !print*,x
            u=matmul(x,dagger(a))/alpha
        case(3) ! SU(3) group
            ! r step
            w=matmul(u,a)
            w1=(real(w(1,1)+w(2,2))*cr+aimag(w(1,1)-w(2,2))*ci)*0.5_r8
            w2=(real(w(1,2)-w(2,1))*cr+aimag(w(1,2)+w(2,1))*ci)*0.5_r8
            detw=dsqrt(dreal(w1*conjg(w1)+w2*conjg(w2)))
            v(1,1)=w1/detw
            v(1,2)=w2/detw
            v(2,1)=-conjg(w2)/detw
            v(2,2)=conjg(w1)/detw
            call hbxmatrix(x,detw)
            y=matmul(x,dagger(v))
            r=ident(nc)
            r(1,1)=y(1,1)
            r(1,2)=y(1,2)
            r(2,1)=y(2,1)
            r(2,2)=y(2,2)
            u=matmul(r,u)

            ! s step
            w=ident(nc)
            w=matmul(u,a)
            w1=(real(w(1,1)+w(3,3))*cr+aimag(w(1,1)-w(3,3))*ci)*0.5_r8
            w2=(real(w(1,3)-w(3,1))*cr+aimag(w(1,3)+w(3,1))*ci)*0.5_r8
            detw=dsqrt(dreal(w1*conjg(w1)+w2*conjg(w2)))
            v(1,1)=w1/detw
            v(1,2)=w2/detw
            v(2,1)=-conjg(w2)/detw
            v(2,2)=conjg(w1)/detw
            call hbxmatrix(x,detw)
            y=matmul(x,dagger(v))
            r=ident(nc)
            r(1,1)=y(1,1)
            r(1,3)=y(1,2)
            r(3,1)=y(2,1)
            r(3,3)=y(2,2)
            u=matmul(r,u)

            ! t step
            w=matmul(u,a)
            w1=(real(w(2,2)+w(3,3))*cr+aimag(w(2,2)-w(3,3))*ci)*0.5_r8
            w2=(real(w(2,3)-w(3,2))*cr+aimag(w(2,3)+w(3,2))*ci)*0.5_r8
            detw=dsqrt(dreal(w1*conjg(w1)+w2*conjg(w2)))
            v(1,1)=w1/detw
            v(1,2)=w2/detw
            v(2,1)=-conjg(w2)/detw
            v(2,2)=conjg(w1)/detw
            call hbxmatrix(x,detw)
            y=matmul(x,dagger(v))
            r=ident(nc)
            r(2,2)=y(1,1)
            r(2,3)=y(1,2)
            r(3,2)=y(2,1)
            r(3,3)=y(2,2)
            u=matmul(r,u)
        case default
            stop
        end select
    endsubroutine sunhb

    ! subroutine that makes one sweep on the lattice via the metropolis algorhitm
    subroutine metropolis(u,s,epsilon,acc)
        complex(kind=r8), intent(inout) :: u(nr,nr,nr,nt,4,nc,nc)
        complex(kind=r8) :: uprime(nc,nc), udif(nc,nc), plaquette(nc,nc), umini(nc,nc)
        real(kind=r8), intent(inout) :: s
        real(kind=r8) :: epsilon, eta, ds
        integer(kind=i4) :: e1, e2, e3, e4, mi, imc
        integer(kind=i4), intent(inout) :: acc
        
        ! we acess each lattice site
        do e1=1,nr
            do e2=1,nr
                do e3=1,nr
                    do e4=1,nt
                        ! on each site we try to evolve each of the four directionr
                        do mi=1,4
                            umini=stample(u,e1,e2,e3,e4,mi)
                            ds=0.0_r8
                            uprime=cz
                            udif=cz
                                    
                            ! calculating the new matrix and the change in the action that it causes
                            uprime=matmul(sunxmatrix(epsilon,nc),u(e1,e2,e3,e4,mi,:,:))
                            call keapinsuN(uprime,nc)
                            udif=uprime-u(e1,e2,e3,e4,mi,:,:)
                            plaquette=matmul(udif,umini)
                            ds=-beta*dreal(ctrace(plaquette,nc))/nc
                            ! now we apply the metropolis criteria
                            call random_number(eta)
                            if((ds.lt.0.0_r8))then
                                u(e1,e2,e3,e4,mi,:,:)=uprime
                                s=s+ds/(6*nt*nr**3.0_r8)
                                acc=acc+1
                            elseif((ds.gt.0.0_r8).and.(exp(-ds).gt.eta))then
                                u(e1,e2,e3,e4,mi,:,:)=uprime
                                s=s+ds/(6*nt*nr**3.0_r8)
                                acc=acc+1
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endsubroutine metropolis

    ! subroutine that makes one sweep on the lattice via the bath algorhitim
    subroutine heatbath(u,nc)
        complex(kind=r8), intent(inout) :: u(nr,nr,nr,nr,4,nc,nc)
        complex(kind=r8) :: stamp(nc,nc) !v(nr,nr,nr,nr,4,nc,nc),
        integer(kind=i4) :: i, j, e1, e2, e3, e4, mi, nc
        do e1=1,nr
            do e2=1,nr
                do e3=1,nr
                    do e4=1,nt
                        do mi=1,4
                            stamp=stample(u,e1,e2,e3,e4,mi) 
                            call sunhb(u(e1,e2,e3,e4,mi,:,:),stamp,nc)
                            call keapinsuN(u(e1,e2,e3,e4,mi,:,:),nc)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endsubroutine heatbath

    ! subroutine that makes one sweep on the lattice via the overrelaxation algorhtim
    subroutine overrelaxation(u,nc)
        complex(kind=r8), intent(inout) :: u(nr,nr,nr,nr,4,nc,nc)
        complex(kind=r8) :: v(nc,nc), unew(nc,nc)
        real(kind=r8) :: ds, zeta
        integer(kind=i4) :: i, j, e1, e2, e3, e4, mi, nc
        do e1=1,nr
            do e2=1,nr
                do e3=1,nr
                    do e4=1,nt
                        do mi=1,4
                            ! we compute the sum of the stamples and project it on the su(N) group
                            v=stample(u,e1,e2,e3,e4,mi)
                            call keapinsuN(v,nc)
                            v=dagger(v)
                            ! now, we try to change the u link
                            unew=matmul(v,dagger(u(e1,e2,e3,e4,mi,:,:)))
                            unew=matmul(unew,v)

                            ! now, we apply a metropolis-like step o acept or not this change
                            ds=-beta*dreal(ctrace(matmul(unew-u(e1,e2,e3,e4,mi,:,:),&
                            &stample(u,e1,e2,e3,e4,mi)),nc))/nc
                            if(ds.le.0.0_r8)then
                                u(e1,e2,e3,e4,mi,:,:)=unew
                            else
                                call random_number(zeta)
                                if(exp(-ds).gt.zeta)then
                                    u(e1,e2,e3,e4,mi,:,:)=unew
                                endif
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endsubroutine overrelaxation
endmodule thermalize