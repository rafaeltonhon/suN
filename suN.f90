program suN
    use lattice
    use thermalize
    use measures
    use mathtools
    implicit none

    ! defining our variables
    real(kind=r8) ::  action
    integer(kind=i4) :: nmc, imix, imc, a, b, i
    
    ! reading the input values and allocating the memory
    open(unit=1,file='suN-in.dat')
    read(1,*) nr, nt    ! lattice size space time
    read(1,*) beta      ! beta
    read(1,*) nmc       ! number of monte carlo sweeps
    read(1,*) a, b
    read(1,*) nc     ! number of colors
    close(1)

    allocate(u(nr,nr,nr,nt,4,nc,nc),w(a,b))

    open(unit=2,file='suN-plaquette.dat')
    call initialize(u,action,1)
    imix=1
    do imc=1,nmc
        ! evolve the system and make the measuments
        if(imix.ne.4) call heatbath(u,nc)
        !if(imix.eq.4)then
        !    call overrelaxation(u,nc)
        !    imix=0
        !endif
        imix=imix+1
        ! smeare the configurations
        !if(i.gt.50) call smearing(u,1.0_r8,200)
        call measurewilson(u,a,b,w)

        ! get away with the data
        write(2,*) imc,w(1,1)
        do i=1,a
            write(100+i,*) imc, w(i,:)
        enddo
    enddo
    close(2)
    deallocate(u)
end program suN