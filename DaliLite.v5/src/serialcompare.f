        program serialcompare
        use soap
        use parsi
        use wolf
        use dalicon
        use filter95
        use dp
        implicit none
        integer outputunit,tmpunit,i1,j2
        character(len=80) dalidatpath_1,dalidatpath_2
        character(len=100000) line
        character(len=16) method
        real rcut,neiborcutoff
        integer maxiter,ndat,dat(64000),nq,nt,n
        logical lfitz,lfirstonly
        character(len=5) cd1,cd2,oldcd1,oldcd2
        character(len=6) key
        character(len=5) qlist(1000000),tlist(1000000)

        if(iargc().lt.3) 
     $          stop "# USAGE: serialcompare DAT_1 DAT_2 method"

        ! input parameters
        call getarg(1,dalidatpath_1)
        call getarg(2,dalidatpath_2)
        call getarg(3,method)
        rcut=4.0
        maxiter=20
        neiborcutoff=12.0
        lfitz=.false.
        lfirstonly=.true.
        outputunit=101
        tmpunit=301
        if(method.eq.'PARSI') call init_parsi()
        if(method.eq.'DALICON') call init_dalicon()
        if(method(1:2).eq.'DP') call dpweights()
        oldcd1='?????'
        oldcd2='?????'
        if(method.eq."DALICON") then
                  if(iargc().lt.4) stop
     $  "# USAGE: mpidali DAT_1 DAT_2 DALICON lfitz"
                  call getarg(4,line)
                  read(line,*) lfitz
        end if
        ! list1 x list2
        if(method.eq.'PARSI'.or.method.eq.'WOLF'.or.
     $          method.eq.'SOAP')then
                call getlist('list1',90,qlist,nq)
                call getlist('list2',90,tlist,nt)
        end if
        select case (method)
          case ("DALICON")
                ! process dalicon_input
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
10              read(*,500,end=99) cd1
                if(cd1(1:3).eq.'END') goto 10
                ! cd2 loop
20              read(*,500,end=99) cd2
500             format(a5)
                if(len_trim(cd2).eq.0) goto 20
                if(cd2(1:3).eq.'END') goto 10
                read(*,*,err=20) ndat
                !write(*,*) '# dalicon input',cd1,cd2,ndat
                if(ndat.lt.1) goto 20
                read(*,*,err=20) (dat(i1),i1=1,4*ndat)
                call dowork_dalicon(cd1,cd2,oldcd1,outputunit,
     $    dalidatpath_1,dalidatpath_2,ndat,dat,lfitz)
                goto 20
99              continue
          case('FILTER95')
                key='refine'
                ! input refine-lines from STDIN
30              read(*,'(A)',end=39) line
                if(line(2:7).ne.key.and.line(1:6).ne.key) goto 30
                n=len_trim(line)
                call dowork_filter95(line,n,1.0,4.0,3,
     $    outputunit,tmpunit,dalidatpath_1,dalidatpath_2,oldcd1,oldcd2)
                goto 30
39              continue
          case("DP")
                key='WOLFIT'
                ! input WOLFITZ-lines from STDIN
130             read(*,'(A)',end=139) line
                if(line(2:7).ne.key.and.line(1:6).ne.key) goto 130
                n=len_trim(line)
                call dowork_dp(line(1:n),n,2.0,oldcd1,oldcd2,
     $    outputunit,dalidatpath_1,dalidatpath_2)
                goto 130
139             continue
          case('PARSI12')
                ! input tuples (cd1,cd2) from STDIN
200             read(*,*,end=219) cd1,cd2
                call dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,lfirstonly)
                goto 200
219             continue
          case('PARSI')
                do i1=1,nq
                  cd1=qlist(i1)
                  do j2=1,nt
                        cd2=tlist(j2)
                        call dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,lfirstonly)
                  end do
                end do
          case('WOLF')
                do i1=1,nq
                  cd1=qlist(i1)
                  do j2=1,nt
                        cd2=tlist(j2)
                        call dowork_wolf(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff)
                  end do
                end do
          case('SOAP')
                 do i1=1,nq
                  cd1=qlist(i1)
                  do j2=1,nt
                        cd2=tlist(j2)
                        call dowork_soap(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter)
                  end do
                end do
        end select

        end program serialcompare
c
c----------------------------------------------------------------------
c

