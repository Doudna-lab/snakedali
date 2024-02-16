c       MPI load balancing to compare list1 x list2
c
c       (c) L. Holm 2013
c
c       timings/WOLF/153lA vs. pdb90    -O3
c       original                        3m 32.641 s
c       N=60               25.616 s         8.971 s
c       N=40               32.298 s        10.103 s
c       N=20               58.280 s        14.352 s
c       N=10            1m 44.762 s        26.096 s
c       N=5             3m 16.134 s        54.101 s
c       N=2            15m  6.771 s     3m 31.059 s
c
c       timings/SOAP/153lA vs. pdb90
c       N=80                            1m 14.957 s
c       N=60            4m 24.023 s     1m 17.190 s
c       N=40            5m 44.843 s     1m 27.715 s
c
c       timings/PARSI/153lA vs. pdb90
c       N=60                               37.499 s
c       N=40                               55.759 s
c       N=20                            1m 40.728 s
c       N=10                            3m 14.979 s
c       N=5                             6m 23.815 s
c       N=2
c
c
c       compile (ignore warnings):
c               module add openmpi-x86_64s
c               mpif90 -O3 mpidali_wolf.f wolf_original.f u3b-8.f subfitzfast.f soap4.f util.f //
c                  ran.f parsi-admin.f parsi-align.f parsi-score.f parsi-stack.f parsi-1.f //
c                  -o mpiwolf
c               mpif90 -O3  gagarubber.f clean.f dimple.f gagadccp.f u3b-8.f lean.f testi.f //
c                  ran.f subfitzfast.f gagagene.f gagahash.f gagatool.f triplet.f util.f //
c                  ssap.f -o dalicon
c
c       execute:
c               mpirun -np 10 ./mpiwolf /data/liisa/DAT/ /data/liisa/DAT/ SOAP
c               mpirun -np 10 ./mpiwolf /data/liisa/DAT/ /data/liisa/DAT/ WOLF
c               mpirun -np 10 ./mpiwolf /data/liisa/DAT/ /data/liisa/DAT/ PARSI
c
c
c       master/slave
c       master_loop:
c               do while workitems in list1 x list2
c                       find idle slave
c                       send workitem (cd1,cd2) to slave
c               send shutdown msg to all slaves
c                       receive acknowledgement
c       slave_loop:
c               report for work
c               do while .true.
c                       receive workitem
c                       if workitem.eq.shutdown:
c                               exit loop
c                       do_work_!
c                       output result to fort.rank
c                       report for work
c       do_work_!
c               setup cd1 if cd1.ne.oldcd1
c               setup cd2 if cd2.ne.oldcd2
c               call compare(cd1,cd2,method)
c       perl wrapper
c               prepare 'list1','list2'
c               remove any existing fort.[10+(1..$NCPU)] output files
c               call "mpirun -np $NCPU mpidali dalidatpath_1 dalidatpath_2 METHOD"
c               collate all fort.[10+(1..$NCPU)] output files
c               ... calculate preliminary zscore, filter bad hits ...
c               ... refine in dalicon ... 
c               ... walk ... refine in dalicon ...
c
c       methods
c               t.b.i. dalicon: input preali!
c               soap4: input cd1,cd2
c               wolf_original: input cd1,cd2
c               parsi: input cd1,cd2
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
        module wolf
        implicit none
        include 'parsizes.for'
        character secstr(maxseg),secstrx(maxseg)
        integer nseg
        real midy(3,maxseg+3),diry(3,maxseg+3)
c
        integer nprot1,nprot2,ns
        character*5 list1(maxprot),list2(maxprot)
c
        integer*4 box(-20:20,-20:20,-20:20),link_a(10000),link_c(10000),
     $          link_b(10000)
        real link_from(3,10000),link_to(3,10000)
        integer lastlink,link_next(10000)
        integer i,j,protcount,bestpair(4)
        real xca(3,maxres),yca(3,maxres),x(3,maxres),y(3,maxres),rms
        integer nx,ny,ali(maxres),k,lali,niter,nsegx
        real midx(3,maxseg+3),dirx(3,maxseg+3),neidist(maxseg,maxseg)
        integer nocc,occ(3,60000)

        contains
c
c------------------------------------------------------------------------------
c
        subroutine dowork_wolf(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff)
        implicit none
        character(len=5) cd1,cd2,oldcd1
        integer outputunit,maxiter
        character*80 dalidatpath_1,dalidatpath_2
        real rcut,neiborcutoff
c
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call initgrid(box,20,20,20)
          nocc=0
          write(*,*) '#setup cd1=',cd1
          call setupprotein(cd1,nsegx,midx,dirx,secstrx,90,xca,nx,
     $          neidist,dalidatpath_1)
          call loadprotein(box,20,20,20,nsegx,midx,dirx,ns,neiborcutoff,
     $          link_a,link_b,link_c,
     $          link_from,link_to,lastlink,link_next,neidist,nocc,occ)
          if(nsegx.le.2) return
        end if
        call setupprotein(cd2,nseg,midy,diry,secstr,90,
     $          yca,ny,neidist,dalidatpath_2)
        if(nseg.le.2) return
        call compare(nseg,midy,diry,secstr,cd2,nseg,secstrx,
     $             link_next,link_from,link_to,link_a,link_b,link_c,
     $             box,20,20,20,neiborcutoff,protcount,bestpair,neidist)
c
c       fitz C(alphas) in best-per-protein trial superimpositions
c
        if(bestpair(1).eq.0.or.bestpair(2).eq.0) return

        call preparex(x,bestpair(3),bestpair(4),0,midx,dirx)
        do j=1,nx
                  do k=1,3
                         x(k,j+3)=xca(k,j)
                  end do
        end do
        call twist(x,3+nx)    ! prepare sets twist frame, transform CA
        do j=1,nx
                  do k=1,3
                        x(k,j)=x(k,j+3)
                  end do
        end do
        call preparex(y,bestpair(1),bestpair(2),0,midy,diry)
        do j=1,ny
                  do k=1,3
                        y(k,j+3)=yca(k,j)
                  end do
        end do
        call twist(y,3+ny)    ! prepare sets twist frame, transform CA
        do j=1,ny
                  do k=1,3
                        y(k,j)=y(k,j+3)
                  end do
        end do
        call fitz(x,nx,y,ny,ali,rcut,maxiter,rms,lali,niter)
        write(outputunit,*) 'WOLFITZ ',cd1,cd2,protcount,
     $            (bestpair(j),j=3,4),(bestpair(j),j=1,2),
     $            rms,lali,niter,nx,(ali(j),j=1,nx)

!510     format('WOLFITZ ',2a5,5i5,f10.1,3i5,<maxres>i4)

        end subroutine dowork_wolf
c
c----------------------------------------------------------------------
c
        end module wolf
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
        module soap

        contains
c
c------------------------------------------------------------------------------
c
        subroutine dowork_soap(cd1,cd2,outputunit,oldcd1,
     $         dalidatpath_1,dalidatpath_2,rcut,maxiter,
     $          ndom1,domns1,domseglist1,nres1,ca1)
        implicit none
        include 'parsizes.for'
        character*5 cd1,cd2
        integer nres1,nres2,outputunit
        integer i,j,ali(maxres+maxres),lali,iter,maxiter,icyc
        integer ali1(maxres)
        integer*2 trace(maxres,maxres),table0(maxres,maxres)
        integer oldlali,bestcyc
        real rms,rcut,u(3,3),t(3),ca1(3,maxres),ca2(3,maxres)
        real r,oldw,totw,x1(3,maxres)
c
        integer nprot1,nprot2,idom1,idom2,iprot1,iprot2
        character*5 list1(maxprot),list2(maxprot),oldcd1,oldcd2
        character*10 string
c
        integer ndom1,ndom2,domns1(maxdom),domns2(maxdom)
        integer domseglist1(maxseg,2,maxdom)
        integer domseglist2(maxseg,2,maxdom)
        character*80 dalidatpath_1,dalidatpath_2
c
        integer domnres1,domres1(maxres),domnres2,domres2(maxres)

c
c       (0) load x,y -domains
c
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call setup(cd1,ndom1,domns1,domseglist1,nres1,ca1,
     $          99,dalidatpath_1)
c         file-not-found handle
          if(ndom1.eq.0) return
        end if
        if(cd2.ne.oldcd2) then
          oldcd2=cd2
          call setup(cd2,ndom2,domns2,domseglist2,nres2,ca2,
     $          99,dalidatpath_2)
c         file-not-found handle
          if(ndom2.eq.0) return
        end if
c
!       type *,'innerloop ',cd1,cd2,nres1,nres2
        oldlali=0
        oldw=0.0
        bestcyc=0
c
c       idom1 * idom2
c
        do idom1=1,ndom1
                domnres1=0
                do i=1,domns1(idom1)
                  do j=domseglist1(i,1,idom1),domseglist1(i,2,idom1)
                        domnres1=domnres1+1
                        domres1(domnres1)=j
                  end do
                end do
                do idom2=1,ndom2
                  domnres2=0
                  do i=1,domns2(idom2)
                    do j=domseglist2(i,1,idom2),domseglist2(i,2,idom2)
                        domnres2=domnres2+1
                        domres2(domnres2)=j
                    end do
                  end do
c       init triangulation
c
        call initequal_fast_domains(table0,nres1,nres2,ali,
     $          domnres1,domres1,domnres2,domres2)
!       write(*,520) (ali(i),i=1,nres1)
        call getut_w(nres1,ali,ca1,ca2,u,t,lali,rms,table0,nres1,nres2)
        call transrotate(ca1,nres1,u,t)
!       write(*,*),'initial lali,rms:',lali,rms,totw
c
c       (2) Gaussian penalty
c
        lali=0
        iter=0
        r=30.0
        rms=99.9
        do while(iter.lt.100.and.r.gt.1.0.and.rms.ge.2.0)
                iter=iter+1
                r=0.90*r
                call fill1_fast(nres1,nres2,ca1,ca2,table0,r)
                call soap_nw_fast(nres1,nres2,table0,ali,trace)
!               write(*,520) (ali(i),i=1,nres1+nres2)
                call soap_getut(nres1+nres2,ali,ca1,ca2,u,t,lali,rms,
     $                  table0,nres1,nres2,totw)
!               write(*,*),'Gaussian -- lali,rms:',lali,rms,iter,totw,r
!               call transrotate(ca1,nres1,u,t)
        end do
        if(totw.gt.oldw) then
                oldw=totw
                oldlali=lali
                do i=1,nres1
                        do j=1,3
                                x1(j,i)=ca1(j,i)
                        end do
                end do
                bestcyc=icyc
        end if
                end do  ! idom2
        end do          ! idom1
c
c       (3) fitz penalty -- best alignment
c
        call fitz(x1,nres1,ca2,nres2,ali1,rcut,maxiter,rms,lali,iter)
        string(1:5)=cd1
        string(6:10)=cd2
        if(string(5:5).eq.' ') string(5:5)='_'
        write(outputunit,*) 'WOLFITZ ',string,nres1,nres2,nint(oldw),
     $    oldlali,bestcyc,rms,lali,iter,nres1,(ali1(j),j=1,nres1)

500     format(a5)
510     format(a80)
520     format(20i4)
!530     format('WOLFITZ ',a10,5i5,f10.1,3i5,<maxres>i4)

999     end subroutine dowork_soap
c
c------------------------------------------------------------------------------
c
        end module soap
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
        module parsi
        implicit none
        include 'parsizes.for'
        logical, parameter :: lpreali=.true.
        logical, parameter :: lnul=.true.
        logical, parameter :: ldb=.true.
        logical, parameter :: lseqtl=.true.
        integer startsize
c
        integer iprot1,iprot2,nprot1,nprot2,ier,na2,nb2,checkx(maxseg)
        character*5 list1(maxprot),list2(maxprot)
        integer segment2(maxres2),nseg2,nres2
        integer dist1sum(maxres1,0:maxres1),minseglen(maxseg)
        character secstr(maxseg),secstr2(maxseg)
        integer nres1,nseg,segmentrange(2,maxseg),na1,nb1,ndom
        integer node_child(2,maxdom),upper(maxseg,maxseg),
     $ lower(maxseg,maxseg),ngap(maxseg),dist2sum(maxres2,0:maxres2)
        integer*2 dist(maxres1*maxres1),dist2(maxres2*maxres2)
        integer ss(maxseg,maxseg),domns(maxdom)
        integer domseglist(maxseg,maxdom)
        logical lfix(maxseg,maxdom),ldom(0:maxdom)
        logical lfix1(maxseg,maxdom),lexdone(maxseg,maxseg)
        integer checkrange(2,maxseg),cut(maxdom),start(maxseg,maxseg)
        character*10 string,str1
        integer segmentrange0(2,maxseg),i

        contains
c
c------------------------------------------------------------------------------
c
        subroutine init_parsi()
        implicit none
        call weights
        call fillscoretable
        end subroutine init_parsi
c
c------------------------------------------------------------------------------
c
        subroutine dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $          dalidatpath_1,dalidatpath_2,lfirstonly)
        implicit none
        include 'parsizes.for'
        character*5 cd1,cd2,oldcd1
        character*80 dalidatpath_1,dalidatpath_2
        logical lfirstonly
        integer outputunit
c
c       structure loop
c
        string(1:5)=cd1
        string(6:10)=cd2
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call getstructure1(nseg,segmentrange,na1,nb1,1,1,ndom,
     $          node_child,dist,nres1,cd1,upper,lower,ss,checkrange,
     $          secstr,lfix,lfix1,cut,domns,domseglist,ier,checkx,
     $          dist1sum,segmentrange0,ldb,dalidatpath_1)
          if(ier.ne.0) return
          call setldom(ldom,ndom,domns,startsize)
c
c         ! be sure to run puutos with the same minseglen parameters !
c
          call setminseglen(secstr,minseglen,nseg,6,8)
          call setngap(ngap,segmentrange,minseglen,10,nseg)
        end if
c
c       sequence loop
c
        call getsequence2(nres2,na2,nb2,dist2,string,dist2sum,
     $          segment2,nseg2,secstr2,na1,nb1,ier,ldb,dalidatpath_2)
        if(ier.ne.0) return
        call initlexdonestart(nseg,ss,start,lexdone)
        call align(ndom,node_child,nres2,lnul,nseg,cut,ldom,domns,
     $          domseglist,upper,lower,segmentrange,dist,nres1,dist2,
     $          dist2sum,lfix,lfix1,start,lexdone,lseqtl,ss,string,
     $          segment2,nseg2,checkrange,secstr,secstr2,checkx,
     $          dist1sum,minseglen,ngap,segmentrange0,ldb,lfirstonly,
     $          lpreali,outputunit)

99      return
        end subroutine dowork_parsi
c
c----------------------------------------------------------------------
c
        end module parsi
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
        program mpidali
        implicit none
        include 'mpif.h'
        integer rank,size,ierr,outputunit
        integer status(MPI_STATUS_SIZE)
        character(len=80) dalidatpath_1,dalidatpath_2
        character(len=5) method
        real rcut,neiborcutoff
        integer maxiter,l

        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size,ierr)
        if(size.lt.2) then
          write(*,*) '# ERROR: at least 2 processes needed'
        else
          call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

          ! input parameters
          call getarg(1,dalidatpath_1)
          call getarg(2,dalidatpath_2)
          call getarg(3,method)
          rcut=4.0
          maxiter=20
          neiborcutoff=12.0

          if(rank.eq.0) then ! master is in control
                ! process list1 x list2
                call master_loop()
                ! stop all workers
                call shutdown(size)
          else if(method(1:4).eq.'SOAP') then
                ! slave waits for requests from master
                outputunit=10+rank
                call slave_loop_soap(rank,outputunit,dalidatpath_1,
     $                  dalidatpath_2,rcut,maxiter)
          else if(method(1:4).eq.'WOLF') then
                ! slave waits for requests from master
                outputunit=10+rank
                call slave_loop_wolf(rank,outputunit,dalidatpath_1,
     $                  dalidatpath_2,rcut,maxiter,neiborcutoff)
          else if(method(1:5).eq.'PARSI') then
                ! slave waits for requests from master
                outputunit=10+rank
                call slave_loop_parsi(rank,outputunit,dalidatpath_1,
     $                  dalidatpath_2,.true.)

          end if
        end if
        call MPI_FINALIZE(ierr)

        end program mpidali
c
c----------------------------------------------------------------------
c
        subroutine master_loop()
        implicit none
        include 'mpif.h'
        include 'parsizes.for'

        integer nprot1,nprot2,i,j
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
                        call send_request_to_slave(cd1,cd2)
                end do
        end do

        end subroutine master_loop
c
c----------------------------------------------------------------------
c
        subroutine slave_loop_soap(myrank,outputunit,dalidatpath_1,
     $  dalidatpath_2,rcut,maxiter)
        use soap
        implicit none
        include 'mpif.h'
        include 'parsizes.for'
        integer myrank,outputunit
        integer status(MPI_STATUS_SIZE),ierr,dat(2),tag,l
        character(len=5) cd1,cd2,oldcd1
        character(len=80) dalidatpath_1,dalidatpath_2
        real rcut
        integer maxiter
        ! soap
        integer ndom1,domns1(maxdom),nres1
        integer domseglist1(maxseg,2,maxdom)
        real ca1(3,maxres)


        oldcd1='?????'
        write(*,*) '#hello from slave',myrank

        ! report for work
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        ! loop waiting for work from master
        do while(.true.)
                ! get work unit (cd1, cd2)
                call MPI_RECV(cd1,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                if(status(MPI_TAG).eq.0) then
                        call acknowledge_shutdown(myrank)
                        return ! exit loop
                end if
                call MPI_RECV(cd2,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                !write(*,*) '#query is',cd1,' ',cd2
                ! process work unit; output goes to file fort.outputunit
                call dowork_soap(cd1,cd2,outputunit,oldcd1,
     $                  dalidatpath_1,dalidatpath_2,rcut,maxiter,
     $                  ndom1,domns1,domseglist1,nres1,ca1)
                ! report for work
                call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     $                  ierr)
        end do

        return
        end subroutine slave_loop_soap
c
c----------------------------------------------------------------------
c
        subroutine slave_loop_wolf(myrank,outputunit,dalidatpath_1,
     $  dalidatpath_2,rcut,maxiter,neiborcutoff)
        use wolf
        implicit none
        include 'mpif.h'
        integer myrank,outputunit
        integer status(MPI_STATUS_SIZE),ierr,dat(2),tag,l
        character(len=5) cd1,cd2,oldcd1
        character(len=80) dalidatpath_1,dalidatpath_2
        real rcut,neiborcutoff
        integer maxiter
        ! soap
        integer ndom1,domns1(maxdom),nres1
        integer domseglist1(maxseg,2,maxdom)
        real ca1(3,maxres)

        oldcd1='?????'
        write(*,*) '#hello from slave',myrank

        ! report for work
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        ! loop waiting for work from master
        do while(.true.)
                ! get work unit (cd1, cd2)
                call MPI_RECV(cd1,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                if(status(MPI_TAG).eq.0) then
                        call acknowledge_shutdown(myrank)
                        return ! exit loop
                end if
                call MPI_RECV(cd2,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                !write(*,*) '#query is',cd1,' ',cd2
                ! process work unit; output goes to file fort.outputunit
                call dowork_wolf(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,rcut,maxiter,neiborcutoff)
                ! report for work
                call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     $                  ierr)
        end do

        return
        end subroutine slave_loop_wolf
c
c----------------------------------------------------------------------
c
        subroutine slave_loop_parsi(myrank,outputunit,dalidatpath_1,
     $  dalidatpath_2,lfirstonly)
        use parsi
        implicit none
        include 'mpif.h'
        integer myrank,outputunit
        integer status(MPI_STATUS_SIZE),ierr,dat(2),tag,l
        character(len=5) cd1,cd2,oldcd1
        character(len=80) dalidatpath_1,dalidatpath_2
        logical lfirstonly

        oldcd1='?????'
        call init_parsi()
        write(*,*) '#hello from slave',myrank

        ! report for work
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        ! loop waiting for work from master
        do while(.true.)
                ! get work unit (cd1, cd2)
                call MPI_RECV(cd1,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                if(status(MPI_TAG).eq.0) then
                        call acknowledge_shutdown(myrank)
                        return ! exit loop
                end if
                call MPI_RECV(cd2,5,MPI_CHARACTER,0,MPI_ANY_TAG,
     $                  MPI_COMM_WORLD,status,ierr)
                !write(*,*) '#query is',cd1,' ',cd2
                ! process work unit; output goes to file fort.outputunit
                call dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $    dalidatpath_1,dalidatpath_2,lfirstonly)
                ! report for work
                call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,
     $                  ierr)
        end do

        return
        end subroutine slave_loop_parsi
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
                write(*,*) '#shutting down',i
                call MPI_SEND('SHUT!',5,MPI_CHARACTER,i,0,
     $                  MPI_COMM_WORLD,ierr)
                ! wait for acknowledgement
                call MPI_RECV(0,0,MPI_INTEGER,i,MPI_ANY_TAG,
     $                          MPI_COMM_WORLD,status,ierr)
                write(*,*) '#shutdown ackn. by',status(MPI_SOURCE)
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
        write(*,*) '#shutdown recvd',myrank
        call MPI_SEND(0,0,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

        return
        end subroutine acknowledge_shutdown
c
c----------------------------------------------------------------------
c
        subroutine send_request_to_slave(cd1,cd2)
        implicit none
        include 'mpif.h'
        integer qprot,qseqlen
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

        return
        end subroutine send_request_to_slave
c
c----------------------------------------------------------------------
c

