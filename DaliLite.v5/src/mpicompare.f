        program mpidali
        implicit none
        include 'mpif.h'
        integer rank,size,ierr,outputunit,tmpunit
        integer status(MPI_STATUS_SIZE)
        character(len=80) line,dalidatpath_1,dalidatpath_2
        character(len=16) method
        real rcut,neiborcutoff
        integer maxiter
        logical lfitz

        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size,ierr)
        if(size.lt.2) stop '# ERROR: at least 2 processes needed'
        if(iargc().lt.3) stop "# USAGE: mpidali DAT_1 DAT_2 method"
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        ! input parameters
        call getarg(1,dalidatpath_1)
        call getarg(2,dalidatpath_2)
        call getarg(3,method)
        rcut=4.0
        maxiter=20
        neiborcutoff=12.0
        lfitz=.false.
        if(method.eq."DALICON") then
                  if(iargc().lt.4) stop
     $  "# USAGE: mpidali DAT_1 DAT_2 DALICON lfitz"
                  call getarg(4,line)
                  read(line,*) lfitz
        end if

        if(rank.eq.0) then ! master is in control
                if(method.eq."DALICON") then
                        ! process dalicon_input
                        call master_loop_dalicon()
                else if(method.eq.'FILTER95') then
                        call master_loop_line('refine')
                else if(method.eq.'DP') then
                        call master_loop_line('WOLFIT')
                else if(method.eq.'PARSI12') then
                        call master_loop_list12()
                else
                        ! process list1 x list2
                        call master_loop_list1_list2()
                end if
                ! stop all workers
                call shutdown(size)
        else
                ! slave waits for requests from master
                outputunit=100+rank
                tmpunit=300+rank
                call slave_loop(rank,tmpunit,outputunit,method,
     $  dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff,.true.,
     $  lfitz)
                close(outputunit)
                write(*,*) '# slave finished',rank
        end if

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)

        end program mpidali
c
c----------------------------------------------------------------------
c
        subroutine master_loop_dalicon()
        implicit none
        include 'mpif.h'
        include 'gagasizes.for'

        integer nprot1,nprot2,i,j,npr,pr(maxres*4)
        character(len=5) cd1,cd2

        ! read from STDIN:
        ! cd1
        !       cd2
        !       npr
        !       pr1(2,1..npr)
        !       pr2(2,1..npr)
        !       ... more cd2 blocks ...
        !       END
        ! ... more cd1 blocks ...
        ! END
        ! cd1 loop
10      read(*,500,end=99) cd1
        if(cd1(1:3).eq.'END') goto 10
        ! cd2 loop
20      read(*,500,end=99) cd2
        if(len_trim(cd2).eq.0) goto 20
        if(cd2(1:3).eq.'END') goto 10
        read(*,*,err=20) npr
        write(*,*) '# dalicon input',cd1,cd2,npr
        if(npr.lt.1) goto 20
        read(*,*,err=20) (pr(i),i=1,4*npr)
        call send_request_to_slave(cd1,cd2,npr,pr)
        goto 20

99      return
500     format(a5)

        end subroutine master_loop_dalicon

c
c----------------------------------------------------------------------
c
        subroutine master_loop_line(key)
        implicit none
        include 'mpif.h'

        integer n,islave,ierr,lkey,status(MPI_STATUS_SIZE)
        character(len=100000) line
        character(len=5) s
        character(len=6) key

        write(*,*) '# start master_loop_line ',key
        s='FILTR'
        ! input refine-lines from STDIN
30      read(*,'(A)',end=39) line
c        write(*,*) '#input key:',line(2:7),' / ',(line(2:7).eq.key.or.
c     $          line(1:6).eq.key)
        if(line(2:7).ne.key.and.line(1:6).ne.key) goto 30
        ! find idle slave
        call MPI_RECV(0,0,MPI_INTEGER,MPI_ANY_SOURCE,0,
     $          MPI_COMM_WORLD,status,ierr)
        islave=status(MPI_SOURCE)
        ! send to slave
        n=len_trim(line)
        ! tag must be 1 for work unit; is 0 for shutdown
        call MPI_SEND(s,5,MPI_CHARACTER,islave,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(n,1,MPI_INTEGER,islave,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(line,n,MPI_CHARACTER,islave,1,MPI_COMM_WORLD,
     $          ierr)
        goto 30
        write(*,*) '# end master_loop_line'
39      return

        end subroutine master_loop_line
c
c----------------------------------------------------------------------
c
        subroutine master_loop_list12()
        implicit none
        include 'mpif.h'
        include 'parsizes.for'

        integer nprot1,nprot2,i,j,dummy(2,1)
        character(len=5) list1(maxprot),list2(maxprot),cd1,cd2

        ! input tuples (cd1,cd2) from STDIN
10      read(*,*,end=19) cd1,cd2
        call send_request_to_slave(cd1,cd2,0,dummy)
        goto 10
19      return

        end subroutine master_loop_list12
c
c----------------------------------------------------------------------
c
        subroutine master_loop_list1_list2()
        implicit none
        include 'mpif.h'
        include 'parsizes.for'

        integer nprot1,nprot2,i,j,dummy(2,1)
        character(len=5) list1(maxprot),list2(maxprot),cd1,cd2

        ! list1 x list2
        nprot1=0
        nprot2=0
        call getlist('list1',90,list1,nprot1)
        call getlist('list2',90,list2,nprot2)
        do i=1,nprot1
                cd1=list1(i)
                do j=1,nprot2
                        cd2=list2(j)
                        call send_request_to_slave(cd1,cd2,0,dummy)
                end do
        end do

        end subroutine master_loop_list1_list2
c
c----------------------------------------------------------------------
c
        subroutine slave_loop(myrank,tmpunit,outputunit,
     $  method,dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff,
     $  lfirstonly,lfitz)
        use soap
        use parsi
        use wolf
        use dalicon
        use filter95
        use dp
        implicit none
        include 'mpif.h'
        integer myrank,tmpunit,outputunit
        character(len=16)  method
        integer status(MPI_STATUS_SIZE)
        character(len=5) cd1,cd2,oldcd1,oldcd2
        character(len=80) dalidatpath_1,dalidatpath_2
        character(len=100000) line
        real rcut,neiborcutoff
        logical lfirstonly,lfitz
        integer maxiter,n,dat(64000)
        integer begintime
        integer, parameter :: timeout=3600 ! 1 hour timeout 

        ! report for work
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        ! loop waiting for work from master
        if(method.eq.'PARSI') call init_parsi()
        if(method.eq.'DALICON') call init_dalicon()
        if(method(1:2).eq.'DP') call dpweights()
        oldcd1='?????'
        oldcd2='?????'
        begintime=MPI_WTIME()
        do while(.true.)       
        !do while(begintime-MPI_WTIME().lt.timeout)
          ! always receive cd1|SHUT!, check for tag
           call MPI_RECV(cd1,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
           if(status(MPI_TAG).eq.0) then
                        call acknowledge_shutdown(myrank)
                        return ! exit loop
           end if
                  ! FILTER95 input is refine-lines
          if(method(1:8).eq."FILTER95".or.method(1:2).eq."DP") then
                call MPI_RECV(n,1,MPI_INTEGER,0,
     $    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(line,n,MPI_CHARACTER,0,
     $    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!                write(*,*) '# slave got bytes',n,line(1:n)
                if(method(1:8).eq."FILTER95") then
                        call dowork_filter95(line,n,1.0,4.0,3,
     $    outputunit,tmpunit,dalidatpath_1,dalidatpath_2,oldcd1,oldcd2)
                else
                        call dowork_dp(line,n,2.0,oldcd1,oldcd2,
     $    outputunit,dalidatpath_1,dalidatpath_2)
                end if
                ! report for work
                call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     $                  ierr)
          else ! input to all others is (cd1,cd2) tuples
                call MPI_RECV(cd2,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                !write(*,*) '# slave ',cd1,' ',cd2,' ',myrank,ierr
                ! process work unit; output goes to file fort.outputunit
                select case (method)
                    case ("PARSI12")
                        call dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,lfirstonly)
                    case ("PARSI")
                        call dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,lfirstonly)
                    case ("SOAP")
                        call dowork_soap(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter)
                    case ("WOLF")
                        call dowork_wolf(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff)
                    case ("DALICON")
                        ! receive prealignment
                        call MPI_RECV(n,1,MPI_INTEGER,0,
     $    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(dat,4*n,MPI_INTEGER,0,
     $    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                        call dowork_dalicon(cd1,cd2,oldcd1,outputunit,
     $    dalidatpath_1,dalidatpath_2,n,dat,lfitz)
                    case default
                        write(*,*) '# method is ',method
                        stop "Unknown method "
                        return ! exit loop
                end select
                ! report for work
                call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     $                  ierr)
          end if
        end do

        return
        end subroutine slave_loop
c
c----------------------------------------------------------------------
c
        subroutine shutdown(size)
        implicit none
        include 'mpif.h'
        integer size
        integer dat(2),i,ierr,status(MPI_STATUS_SIZE)

        ! send DIETAG to all workers; wait for acknowledgement

        write(*,*) '# SHUTDOWN!',size
        dat(1)=0
        dat(2)=0
        ! shutdown slaves
        do i=1,size-1
                !write(*,*) '#shutting down',i
                call MPI_SEND('SHUT!',5,MPI_CHARACTER,i,0,
     $                  MPI_COMM_WORLD,ierr)
                ! wait for acknowledgement
                call MPI_RECV(0,0,MPI_INTEGER,i,MPI_ANY_TAG,
     $                          MPI_COMM_WORLD,status,ierr)
                !write(*,*) '#shutdown ackn. by',status(MPI_SOURCE)
        end do

        return
        end subroutine shutdown
c
c----------------------------------------------------------------------
c
        subroutine acknowledge_shutdown(myrank)
        implicit none
        include 'mpif.h'
        integer ierr,myrank

        ! acknowledge shutdown msg
        !write(*,*) '#shutdown recvd',myrank
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        return
        end subroutine acknowledge_shutdown
c
c----------------------------------------------------------------------
c
        subroutine send_request_to_slave(cd1,cd2,npr,pr)
        implicit none
        include 'mpif.h'
        integer npr,pr(4*npr)
        character(len=5) cd1,cd2
        integer islave,status(MPI_STATUS_SIZE),ierr

        ! find idle slave
        call MPI_RECV(0,0,MPI_INTEGER,MPI_ANY_SOURCE,0,
     $          MPI_COMM_WORLD,status,ierr)
        islave=status(MPI_SOURCE)
        !write(*,*) '#slave',islave,' reported for work ',cd1,' ',cd2
        ! send request to slave
        call MPI_SEND(cd1,5,MPI_CHARACTER,islave,1,
     $          MPI_COMM_WORLD,ierr)
        call MPI_SEND(cd2,5,MPI_CHARACTER,islave,1,
     $          MPI_COMM_WORLD,ierr)
        !write(*,*) '# work unit ',cd1,' ',cd2,' ',islave,ierr
        if(npr.eq.0) return ! list1 x list2 methods
        ! DALICON, sends alignment
        call MPI_SEND(npr,1,MPI_INTEGER,islave,1,
     $          MPI_COMM_WORLD,ierr)
        call MPI_SEND(pr,4*npr,MPI_INTEGER,islave,1,
     $          MPI_COMM_WORLD,ierr)

        return
        end subroutine send_request_to_slave
c
c----------------------------------------------------------------------
c

