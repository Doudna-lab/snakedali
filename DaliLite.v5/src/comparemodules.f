c       MPI load balancing to compare list1 x list2
c
c       (c) L. Holm 2013
c
c       timings/WOLF/153lA vs. pdb90    -O3              blask/reverse
c       original                        3m 32.641 s
c       N=60               25.616 s         8.971 s
c       N=40               32.298 s        10.103 s     0m12.854s/ 0m11.357s
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
c       N=60                               37.499 s     0m51.354s
c       N=40                               55.759 s      1m16.615s/0m48.006s
c       N=20                            1m 40.728 s     2m18.581s
c       N=10                            3m 14.979 s     4m31.089s
c       N=5                             6m 23.815 s
c       N=2                            28m 41.915 s
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
c               mpif90 -O1 dp.f util.f -o dp ! crashes with higher -O level !
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
c
c
c       methods
c               t.b.i. dalicon: input preali!
c               soap4: input cd1,cd2
c               wolf_original: input cd1,cd2
c               parsi: input cd1,cd2
c
c>> search Protein Data Bank, pairwise, one-against-all, all-against-all
c>> map PDB-id (or sequence!) to precomputed
c----------------------------------------------------------------------
c----------------------------------------------------------------------
        subroutine compressblocks(nblock,l1,r1,l2,r2)
        implicit none
        include 'parsizes.for'
        integer nblock,l1(maxres),r1(maxres),l2(maxres),r2(maxres)
        integer i,n

        i=0
        n=0
        do while(i.lt.nblock)
                i=i+1
                n=n+1
                !write(*,*) '#compress',i,l1(i),r1(i),l2(i),r2(i)
                l1(n)=l1(i)
                l2(n)=l2(i)
                r1(n)=r1(i)
                r2(n)=r2(i)
                if(i.lt.nblock) then
                  do while(l1(i).eq.r1(i).and.l2(i).eq.r2(i).and.
     $  l1(i+1).eq.r1(i)+1.and.l2(i+1).eq.r2(i)+1)
                        if(i+1.ge.nblock) exit ! loop
                        i=i+1
                  end do
                  r1(n)=r1(i)
                  r2(n)=r2(i)
                end if
                !write(*,*) '#compressed',n,l1(n),r1(n),l2(n),r2(n)
        end do
        nblock=n

        end subroutine compressblocks
c
c------------------------------------------------------------------------------
c
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
        integer*4 box(-20:20,-20:20,-20:20),link_a(100000),
     $          link_b(100000),link_c(100000)
        real link_from(3,100000),link_to(3,100000)
        integer lastlink,link_next(100000)
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
        integer i,nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
c

        if(cd1.ne.oldcd1) then
          ! output daliocon format
!          write(outputunit,500) 'END  '
!          write(outputunit,500) cd1
          oldcd1=cd1
          call initgrid(box,20,20,20)
          nocc=0
c          write(*,*) '#setup cd1=',cd1
          call setupprotein(cd1,nsegx,midx,dirx,secstrx,outputunit+100,
     $          xca,nx,neidist,dalidatpath_1)
c          write(*,*) 'cd1',cd1,nsegx,nx
          call loadprotein(box,20,20,20,nsegx,midx,dirx,ns,neiborcutoff,
     $          link_a,link_b,link_c,
     $          link_from,link_to,lastlink,link_next,neidist,nocc,occ)
          if(nsegx.le.2) return
        end if
c        write(*,*) '#setup cd2=',cd2
        call setupprotein(cd2,nseg,midy,diry,secstr,outputunit+100,
     $          yca,ny,neidist,dalidatpath_2)
c        write(*,*) 'cd2',cd2,nseg,ny
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
        !write(*,*) '#calling fitz1 ',cd1,nx,' ',cd2,ny
        call fitz(x,nx,y,ny,ali,rcut,maxiter,rms,lali,niter)
        ! conert ali to blocks
        nblock=0
        do i=1,nx
                  if(ali(i).ne.0) then
                        nblock=nblock+1
                        l1(nblock)=i
                        r1(nblock)=i
                        l2(nblock)=ali(i)
                        r2(nblock)=ali(i)
                  end if
        end do
        ! output dalicon format
c        write(outputunit,500) cd2,'*'
c        write(outputunit,*) nblock
!        write(*,*) '#l1wolf',(l1(i),r1(i),i=1,nblock)
!        write(*,*) '#l2wolf',(l2(i),r2(i),i=1,nblock)
c       write(outputunit,510) cd1,cd2,1,0,nblock,
        call compressblocks(nblock,l1,r1,l2,r2)
        write(outputunit,*) 'WOLFITZ ',cd1,cd2,nblock,
     $          (l1(i),r1(i),i=1,nblock),
     $          (l2(i),r2(i),i=1,nblock)
!        write(outputunit,*) 'WOLFITZ ',cd1,cd2,protcount,
!     $            (bestpair(j),j=3,4),(bestpair(j),j=1,2),
!     $            rms,lali,niter,nx,(ali(j),j=1,nx)

!510     format('WOLFITZ ',2a5,5i5,f10.1,3i5,<maxres>i4)
500     format(a5,a1)
510     format('refine',2a5,i4,i20,i4,800i4)

        end subroutine dowork_wolf
c
c----------------------------------------------------------------------
c
        end module wolf
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
        module filter95
        implicit none
        include 'parsizes.for'
        integer nseg,ali(4*maxres),ndom
        integer node_size(maxdom)
        integer minscore(2,maxdom),nres,nres2,ndom2,node_size2(maxdom)
        real ca(3,maxres),ca2(3,maxres)
        logical lkeep(maxdom)
        character node_type(maxdom),node_type2(maxdom)
        integer*2 d1(maxres,maxres),d2(maxres,maxres)
        integer node_nseg(maxdom),node_nseg2(maxdom)
        integer segmentrange(2,maxseg,maxdom)
        integer segmentrange2(2,maxseg,maxdom)
 
        private
        public dowork_filter95
 
        contains
c
c----------------------------------------------------------------------------
c
        subroutine dowork_filter95(line,n,zcut1,fitzrcut,fitzmaxiter,
     $          outunit,tmpunit,dalidatpath_1,dalidatpath_2,
     $          oldcd1,oldcd2)
        implicit none
        include 'parsizes.for'
        real zcut1
        integer outunit,tmpunit
        character(len=100000) line
c
        character*5 oldcd1,oldcd2,cd1,cd2
        integer idom,score,n,i,j,k,l
        real fitzrcut
        integer fitzmaxiter
        character*80 constructfilnam,filnam,dalidatpath_1,dalidatpath_2
        real mean,sigma,x,rms,u(3,3),t(3)
        integer lali,niter,fitzali(maxres),xscore,resix(maxres)
        integer xiser(maxres)
        integer nx,oldidom
        real xca(3,maxres),zscore
        logical lnew1,ldebug
        parameter(ldebug=.false.)
        integer tmpali(maxres)
c
        oldidom=0
10      continue
!        write(*,*) '#this is dowork_filter95',line(1:n)
        read(line(8:12),*,end=19) cd1
        read(line(13:17),*,end=19) cd2
        read(line(18:n),*,end=19) idom,score,nseg,(ali(i),i=1,nseg*4)
!        write(*,*),'input:', oldcd1,' ',cd1,' ',oldcd2,' ',cd2,
!     $          idom,lkeep(idom),score,nseg,(ali(i),i=1,nseg*4)

        if(cd1.eq.oldcd1.and.cd2.eq.oldcd2.and.idom.eq.oldidom) return
        lnew1=.false.
        if(cd1.ne.oldcd1) then
                filnam=constructfilnam(cd1,dalidatpath_1,'.dat')
                open(tmpunit,file=filnam,status='old',err=11)
                goto 12
                ! try second dalidatpath
11              filnam=constructfilnam(cd1,dalidatpath_2,'.dat')
                open(tmpunit,file=filnam,status='old',err=10)
12              continue
                call readproteindata95(ndom,node_type,node_size,nres,
     $                  ca,segmentrange,node_nseg,tmpunit)
                close(tmpunit)
                call getdist95(ca,min(maxres,nres),d1)
                do i=1,ndom
                  lkeep(i)=(node_type(i).eq.'*'.or.node_type(i).eq.'+')
                  if(ldebug) write(*,*) 'idom,lkeep ',i,lkeep(i)
                end do
                !
                ! minscore(1,2) is Z-score mean,sigma using query domain
                ! size
                !
                do i=1,ndom
                  x=min(float(node_size(i)),400.0)
                  mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
                  sigma=max(1.0,0.50*mean)
                  minscore(1,i)=nint(10000*mean)
                  minscore(2,i)=nint(10000*sigma)
                end do
                oldcd1=cd1
                lnew1=.true.
        end if
        if(.not.lkeep(idom)) return
!        ! keep directly if Z-score > 2.0
!        zscore=float(score-minscore(1,idom))/minscore(2,idom)
!        if(zscore.ge.zcut1) then
!          write(outunit,*) zscore,cd1,cd2,idom,score,nseg,
!     $          (ali(i),i=1,nseg*4)
!          if(ldebug)write(*,*),'save 1:',idom,score,zscore,' ',cd1,cd2
!          oldcd2=cd2
!          oldidom=idom
!          return
!        end if
        !
        ! give lower scoring pairs a chance using fitz
        !
        if(idom.ne.oldidom.or.lnew1) then
                do i=1,min(maxres,nres)
                        resix(i)=0
                        xiser(i)=0
                end do
                nx=0
                do i=1,node_nseg(idom)
                if(ldebug) write(*,*)
     $                  'domain:',idom,i,(segmentrange(j,i,idom),j=1,2)
                   do j=segmentrange(1,i,idom),segmentrange(2,i,idom)
                        if(j.gt.0.and.j.le.maxres) then
                                nx=nx+1
                                resix(nx)=j
                                xiser(j)=nx
                                do k=1,3
                                        xca(k,nx)=ca(k,j)
                                end do
                        end if
                   end do
                end do
                oldidom=idom
                if(ldebug.and.nx.ne.node_size(idom)) write(*,*),
     $                  'ERROR: nx.ne.node_size !',cd1,cd2,
     $                  idom,nx,node_size(idom),node_nseg(idom),
     $                  (segmentrange(1,i,idom),segmentrange(2,i,idom),
     $                  i=1,node_nseg(idom))
        end if
!        write(*,*) '#checkpoint0 ',cd2,' ',oldcd2,nres2
        if(cd2.ne.oldcd2) then
                filnam=constructfilnam(cd2,dalidatpath_2,'.dat')
                open(tmpunit,file=filnam,status='old',err=10)
                call readproteindata95(ndom2,node_type2,node_size2,
     $            nres2,ca2,segmentrange2,node_nseg2,tmpunit)
                close(tmpunit)
!                write(*,*) '#after readproteindata95 ',cd2,nres2
                call getdist95(ca2,nres2,d2)
!                write(*,*) '#after getdist95 ',cd2,nres2
                oldcd2=cd2
                !write(*,*),cd1,cd2
        end if
!        write(*,*) '#checkponit1',nres2
        if(ldebug)write(*,*),'ali',(ali(i),i=1,nseg*4)
        do i=1,nx
                fitzali(i)=0
        end do
!        write(*,*) '#before l',nres2
        do i=1,nseg
                k=(i-1)*2
                if(ali(k+1).gt.0) then
                        do j=ali(k+1),ali(k+2)
                                l=xiser(j)
                                if(l.gt.0.and.l.le.nres) then
                                        fitzali(l)=
     $                                    ali(k+nseg*2+1)+j-ali(k+1)
!        write(*,*) '#l',l,k,j,ali(k+nseg*2+1),ali(k+1),fitzali(l),nres2
                                end if
                        end do
                end if
        end do
!        write(*,*),'fitzali1',(fitzali(i),i=1,nx),cd1,nx,' ',cd2,nres2
        call getut(nx,fitzali,xca,ca2,u,t,lali,rms)
!        write(*,*),'getut done',lali,rms,cd1,nx,' ',cd1,nx,' ',cd2,nres2
        call transrotate(xca,nx,u,t)
!        write(*,*),'transrotate done ',cd1,nx,' ',cd2,nres2
!        write(*,*) '#calling fitz2 ',cd1,nx,' ',cd2,nres2
        call fitz(xca,nx,ca2,nres2,fitzali,fitzrcut,fitzmaxiter,rms,
     $          lali,niter)
        if(ldebug)write(*,*),'fitzali2',(fitzali(i),i=1,nx)
        do i=1,min(maxres,nres)
                tmpali(i)=0
        end do
        do i=1,nx
                if(fitzali(i).gt.0) then
                        tmpali(resix(i))=fitzali(i)
                end if
        end do
        call gettotscore95(tmpali,d1,nres,d2,nres2,x)
        if(ldebug)write(*,*),'gettotscore95 done',x
        xscore=nint(10000*x)
        if(ldebug)write(*,*),'fitz done',rms,lali,niter,xscore,score
        if(xscore.gt.score) then        ! refine fitzed ali
                nseg=0
                do i=1,nx
                        if(fitzali(i).gt.0) nseg=nseg+1
                end do
                j=0
                do i=1,nx
                        if(fitzali(i).gt.0) then
                                j=j+1
                                ali((j-1)*2+1)=resix(i)
                                ali((j-1)*2+2)=resix(i)
                                ali((j-1)*2+nseg*2+1)=fitzali(i)
                                ali((j-1)*2+nseg*2+2)=fitzali(i)
                        end if
                end do
                score=xscore
        end if
        zscore=float(score-minscore(1,idom))/minscore(2,idom)
        if(ldebug) write(*,*),'#zscore',idom,score,zscore,' ',cd1,cd2
        if(zscore.lt.zcut1) return ! reject low-scoring even after fitz
        write(outunit,*) zscore,cd1,cd2,idom,score,nseg,
     $          (ali(i),i=1,nseg*4)
19      return 

500     format(6x,2a5,i4,i20,i4,16000i4)
510     format(f5.1,1x,2a5,i4,i20,i4,16000i4)
520     format(a5,a1,a1)

        return
        end subroutine dowork_filter95
c
c----------------------------------------------------------------------------
c
        subroutine getdist95(ca,nres,d1)
        implicit none
        include 'parsizes.for'
        integer nres
        real ca(3,maxres)
        integer*2 d1(nres,nres)
c
        integer i,j
c
        do i=1,min(maxres,nres)
                d1(i,i)=0.0
                do j=i+1,min(maxres,nres)
                  d1(i,j)=intdistance(ca(1,i),ca(2,i),ca(3,i),
     $                   ca(1,j),ca(2,j),ca(3,j))
                  d1(j,i)=d1(i,j)
                end do
        end do

        return
        end subroutine getdist95
c
c----------------------------------------------------------------------------
c
        subroutine readproteindata95(ndom,node_type,node_size,nres,ca,
     $          segmentrange,node_nseg,iunit)
        implicit none
        include 'parsizes.for'
        integer ndom,nres,nseg,iunit,segmentrange(2,maxseg,maxdom)
        integer node_size(maxdom) ,node_nseg(maxdom)
        character node_type(maxdom)
        real ca(3,maxres),x
c
        integer i,j,k,idom,iseg,n
c
        read(iunit,500) nres,nseg
!        write(*,*) '#readproteindata95',nres,nseg,maxres,maxseg
        do iseg=1,nseg
                read(iunit,510) i
        end do
        n=nres-maxres
        if(nres.gt.maxres) nres=maxres
        read(iunit,520) ((ca(j,i),j=1,3),i=1,nres),(x,i=1,3*n)
        read(iunit,500) ndom
        do idom=1,ndom
                read(iunit,530) i
        end do
        read(iunit,500) ndom
        do idom=1,min(maxdom,ndom)
          read(iunit,530) i,node_type(i),j,j,node_size(i),node_nseg(i),
     $    ((segmentrange(j,k,i),j=1,2),k=1,node_nseg(i))
        end do

500     format(10x,4i5,2x,64000a1)
510     format(6i10)
520     format(10f8.1)
530     format(i4,1x,a1,1x,4i4,64000i4)

        return
        end subroutine readproteindata95
c
c-------------------------------------------------------------------------------
c
        subroutine gettotscore95(ali1,d1,nres1,d2,nres2,totscore)
        implicit none
        include 'parsizes.for'
        real totscore
        integer nres1,nres2
        integer ali1(maxres)
c
        integer*2 d1(nres1,nres1),d2(nres2,nres2)
c
        integer i,j,k,l,q,r,a(maxres),n
        real x
c
        totscore=0.0
        n=0
        do i=1,min(maxres,nres1)
                if(ali1(i).gt.0) then
                        n=n+1
                        a(n)=i
                end if
        end do
        do i=1,n
                k=a(i)
                do j=1,n
                        l=a(j)
                        q=abs(ali1(k))
                        r=abs(ali1(l))
                        if(q.lt.1.or.r.lt.1) cycle ! HACK 
                        if(q.gt.nres2.or.r.gt.nres2) cycle ! HACK
                        x=scorefun95(d1(k,l),d2(q,r))
                        totscore=totscore+x
                end do
        end do
c
500     format(20i4)
c
        return
        end subroutine gettotscore95
c
c-------------------------------------------------------------------------------
c
        real function scorefun95(r1,r2) result(s)
        implicit none
        integer*2 r1,r2
        real r
c
        r=float(r1+r2)/200.0
        s=(r1-r2)/100.0
        if(r.gt.0.01) then
                s=(0.20-abs(s)/r)*exp(-r*r/400.0)
        else
                s=0.20
        end if
c
        return
        end function scorefun95
c
c-------------------------------------------------------------------------------
c
        integer*2 function intdistance(v1,v2,v3,u1,u2,u3) result(d)
        implicit none
        real v1,v2,v3,u1,u2,u3,x1,x2,x3
c
        x1=v1-u1
        x2=v2-u2
        x3=v3-u3
        d=nint(100.0*sqrt(x1*x1+x2*x2+x3*x3))
c
        return
        end function intdistance
c
c----------------------------------------------------------------------
c

        function distance(a1,a2,a3,b1,b2,b3)
        implicit none
        real distance,a1,a2,a3,b1,b2,b3
c
        distance=sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))
c
        return
        end function distance
c
c------------------------------------------------------------------------------
c


        end module filter95
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
        module soap
        include 'parsizes.for'
        integer ndom1,domns1(maxdom),nres1
        integer domseglist1(maxseg,2,maxdom)
        real ca1(3,maxres),ca2(3,maxres)
        integer ndom2,domns2(maxdom),nres2
        integer domseglist2(maxseg,2,maxdom)

        contains
c
c------------------------------------------------------------------------------
c
        subroutine dowork_soap(cd1,cd2,outputunit,oldcd1,
     $         dalidatpath_1,dalidatpath_2,rcut,maxiter)
        implicit none
        character*5 cd1,cd2
        integer outputunit
        integer i,j,ali(maxres+maxres),lali,iter,maxiter,icyc
        integer ali1(maxres)
        integer*2 trace(maxres,maxres),table0(maxres,maxres)
        integer oldlali,bestcyc
        real rms,rcut,u(3,3),t(3)
        real r,oldw,totw,x1(3,maxres)
c
        integer nprot1,nprot2,idom1,idom2,iprot1,iprot2
        character*5 list1(maxprot),list2(maxprot),oldcd1,oldcd2
        character*10 string
c
        character*80 dalidatpath_1,dalidatpath_2
c
        integer domnres1,domres1(maxres),domnres2,domres2(maxres)
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)

c
c       (0) load x,y -domains
c
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call setup(cd1,ndom1,domns1,domseglist1,nres1,ca1,
     $          99,dalidatpath_1,dalidatpath_2)
c         file-not-found handle
          if(ndom1.eq.0) return
        end if
        if(cd2.ne.oldcd2) then
          oldcd2=cd2
          call setup(cd2,ndom2,domns2,domseglist2,nres2,ca2,
     $          99,dalidatpath_1,dalidatpath_2)
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
                do i=1,min(maxres,nres1)
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
        !write(*,*) '#calling fitz3 ',cd1,nres1,' ',cd2,nres2
        call fitz(x1,nres1,ca2,nres2,ali1,rcut,maxiter,rms,lali,iter)
        string(1:5)=cd1
        string(6:10)=cd2
        if(string(5:5).eq.' ') string(5:5)='_'
        write(*,*) '#soap ali1:',cd1,cd2,rms,lali,nres1,nres2,
     $          (ali1(i),i=1,nres1)
        nblock=0
        do i=1,min(maxres,nres1)
                  if(ali1(i).ne.0) then
                        nblock=nblock+1
                        l1(nblock)=i
                        r1(nblock)=i
                        l2(nblock)=ali1(i)
                        r2(nblock)=ali1(i)
                  end if
        end do
        call compressblocks(nblock,l1,r1,l2,r2)
        write(outputunit,*) 'WOLFITZ ',cd1,cd2,nblock,
     $          (l1(i),r1(i),i=1,nblock),
     $          (l2(i),r2(i),i=1,nblock)
c        write(outputunit,*) 'WOLFITZ ',string,nres1,nres2,nint(oldw),
c     $    oldlali,bestcyc,rms,lali,iter,nres1,(ali1(j),j=1,nres1)

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
        logical, parameter :: lpreali=.false.
        logical, parameter :: lnul=.true.
        logical, parameter :: ldb=.true. ! read inputs from list
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
        call weights()
        call fillscoretable()
        end subroutine init_parsi
c
c------------------------------------------------------------------------------
c
        subroutine dowork_parsi(cd1,cd2,outputunit,oldcd1,
     $          dalidatpath_1,dalidatpath_2,lfirstonly)
        implicit none
ccc        include 'parsizes.for'
        character*5 cd1,cd2,oldcd1
        character*80 dalidatpath_1,dalidatpath_2
        logical lfirstonly
        integer outputunit
c
c       structure loop
c
        string(1:5)=cd1
        string(6:10)=cd2
        lskipflag=.false. ! common
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call getstructure1(nseg,segmentrange,na1,nb1,1,1,ndom,
     $          node_child,dist,nres1,cd1,upper,lower,ss,checkrange,
     $          secstr,lfix,lfix1,cut,domns,domseglist,ier,checkx,
     $          dist1sum,segmentrange0,ldb,dalidatpath_1)
          write(*,*) '# getstructure1 ',cd1,' ',cd2,' ',ier
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
        write(*,*) '# getsequence2 ',cd1,' ',cd2,' ',ier
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
        module dalicon
        implicit none
        include 'gagasizes.for'
        character chainid1,chainid2,seq1(maxres0),seq2(maxres0),
     $          struc1(maxres0),struc2(maxres0)
        character*4 code1,code2,resno1(maxres0),resno2(maxres0)
        integer nres1,nres2,relacc1(maxres0),relacc2(maxres0)
        integer*2 d1(maxres,maxres),d2(maxres,maxres)
        real ca1(3,maxres0),ca2(3,maxres0)
        character*80 hdr1,hdr2,infile
c
        integer ierr,outunit
        logical file_exists
c
        integer*2 genepool(5,maxpair)
        integer ngene
        integer*2 ali1(maxres,popsize),preali1(maxres0)
        real score(popsize),oldscore,ds
        integer i,j,k,npr,pr1(2,maxres),pr2(2,maxres)
        logical ldef,lpreali,ltop0
        integer ntetra,nres0,tsil(maxres),nres02,tsil2(maxres)
        integer*2 tetrapool(2,maxpair)
c
        integer ali(maxres)
        real u(3,3),t(3),fitzrcut,rms
        integer nx,fitzmaxiter,nali,niter
        parameter (fitzrcut=4.0, fitzmaxiter=10)
c
        integer jnk1(dim1)
        common /junk1/jnk1
        integer jnk23(dim23)
        common /junk23/jnk23
        integer jnk4(dim4)
        common /junk4/jnk4
        integer jnk5(dim5)
        common /junk5/jnk5
        integer jnk6(dim6)
        common /junk6/jnk6


        contains
c
c----------------------------------------------------------------------
c
        subroutine init_dalicon()
        implicit none

        nres0=maxres
        ldef=.true.
        open(90,file='dali.default',status='old',err=5)
        goto 6
5       ldef=.false.
        write(*,*) ' *** dali.default not found *** that''s ok ***'
6       call defaults(90,ldef)
        write(*,*),'defaults done'
        close(90)
c
        call gagaweights()
        do i=1,maxres
                tsil(i)=0
        end do

        end subroutine init_dalicon
c
c----------------------------------------------------------------------
c
        subroutine dowork_dalicon(cd1,cd2,oldcd1,outputunit,
     $          dalidatpath_1,dalidatpath_2,npr,pr,lfitz)
        implicit none
        character(len=5) cd1,cd2,oldcd1
        integer outputunit,npr,pr(4*npr),i,j,offset
        character(len=80) dalidatpath_1,dalidatpath_2
        logical lfitz
        integer iter
        integer, parameter :: maxiter=100
        logical, parameter :: ldebug=.false.

        ! setup cd1
        if(cd1.ne.oldcd1) then
                call setup_new(cd1,nres1,ca1,99,dalidatpath_1,seq1)
                call getgagadist(ca1,nres1,d1,maxres)
                oldcd1=cd1
        end if
        do i=1,min(maxres,nres1)
                tsil(i)=i
                preali1(i)=0
        end do
        ! setup cd2
        call setup_new(cd2,nres2,ca2,99,dalidatpath_2,seq2)
        call getgagadist(ca2,nres2,d2,maxres)
        do j=1,min(maxres,nres2)
                tsil2(j)=j
        end do
        code1=cd1(1:4)
        chainid1=cd1(5:5)
        code2=cd2(1:4)
        chainid2=cd2(5:5)
c
c       convert prealignment ranges pr1,pr2 to preali1
c
        nx=0
        do i=1,npr,2
                k=pr(npr*2+i)
                do j=pr(i),pr(i+1)
                        preali1(j)=k
                        k=k+1
                        nx=nx+1
                end do
        end do
        if(ldebug) write(*,*),'preali1 done',k,(tsil(j),j=1,nres1),
     $          (tsil2(j),j=1,nres2)
c
c       run fitz to extend short prealignment
c
        if(lfitz) then
!               write(*,*) 'initial preali1',nres1,(preali1(i),i=1,nres1)
c               preali1 is integer*2, getut expects ali as integer
                do i=1,min(maxres,nres1)
                        ali(i)=preali1(i)
                end do
c               initial superimposition by u3b
                call getut(nres1,ali,ca1,ca2,u,t,nali,rms)
!               write(*,*) 'getut result:',nali,rms,nres1,(ali(i),i=1,nres1)
                call transrotate(ca1,nres1,u,t)
c               fitz iterations
                if(ldebug)write(*,*) '#fitz4 ',cd1,nres1,' ',cd2,nres2
                call fitz(ca1,nres1,ca2,nres2,ali,fitzrcut,fitzmaxiter,
     $                  rms,nali,niter)
c               overwrite preali1 if fitzed is longer
                if(nali.gt.nx) then
                        do i=1,min(maxres,nres1)
                                preali1(i)=ali(i)
!                       write(*,*) 'fitz result ',i,ali(i),preali1(i)
                        end do
                        if(ldebug) write(*,*) 'preali1',rms,nali,niter,
     $                          (preali1(i),i=1,nres1)
                end if
        end if
c
c       reverse Polish stack mode of operation
c
            call compresspreali(preali1,tsil,tsil2,nres1,nres2,nres0)
            if(ldebug) write(*,*),'compresspreali done'
            call testi(nres1,nres2,preali1,d1,d2,ntetra,tetrapool)
            if(ldebug) write(*,*),'testi done',nres1,nres2,ntetra
            lverb=.false.
            score(1)=0.0
            do i=1,maxres
                ali1(i,1)=0
            end do
            if(ntetra.gt.0) call lean_mc(50,4,nres1,nres2,d1,d2,.true.,
     $          preali1,itrim(2),score,ali1,ntetra,tetrapool)
            write(*,*),'lean_mc done', score(1)
c
c           iterate until improvement < 10
c
            oldscore=score(1)
            ds=oldscore
            iter=0
            do while(ds.gt.refitol.and.iter.lt.maxiter)
              iter=iter+1
              do i=1,maxres
                preali1(i)=ali1(i,1)
              end do
              call testi(nres1,nres2,preali1,d1,d2,ntetra,tetrapool)
              lverb=.false.
              score(1)=0.0
              do i=1,maxres
                ali1(i,1)=0
              end do
              if(ntetra.gt.0) call lean_mc(50,4,nres1,nres2,d1,d2,
     $          .true.,preali1,itrim(2),score,ali1,ntetra,tetrapool)
              ds=score(1)-oldscore
              oldscore=score(1)
              if(ldebug) write(*,*) '#iter',iter,ds,score(1)
            end do
            if(ldebug) write(*,*) (ali1(k,1),k=1,nres1)
c
            if(ldebug) write(*,*) 'call output',(score(i),i=1,popsize)
            call output(code1,chainid1,code2,chainid2,nres1,nres2,seq1,
     $  seq1,struc1,struc2,resno1,resno2,1,hdr1,hdr2,ca1,ca2,
     $  relacc1,relacc2,d1,d2,ali1,score,nres0,tsil,nres02,tsil2,
     $  outputunit)

        return

500     format(/' comparing: ',a4,a1,i5,' res to ',a4,a1,i5,' res ',l5)
510     format(a20,' output file exists already !!!')
520     format(' antiparallel matches are ',a20/)

        end subroutine dowork_dalicon
c
c----------------------------------------------------------------------
c
        end module dalicon
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
        module dp
        implicit none
        include 'parsizes.for'
        integer i,ndom1,ndom2,domns1(maxdom),domns2(maxdom),nres1,nres2
        integer domseglist1(maxseg,2,maxdom)
        integer domseglist2(maxseg,2,maxdom)
        real ca1(3,maxres),ca2(3,maxres),wght(0:100)
        integer*2 d1(maxres,maxres),d2(maxres,maxres)
        character seq1(maxres),seq2(maxres)

        private
        public dowork_dp,dpweights


        contains

c------------------------------------------------------------------------------
c
        subroutine dowork_dp(line,n,zcut,oldcd1,oldcd2,outputunit,
     $          dalidatpath_1,dalidatpath_2)
        implicit none
        include 'parsizes.for'
        character*5 cd1,cd2,oldcd1,oldcd2
        character(len=10) cd1cd2
        integer outputunit,i,n
        character*80 dalidatpath_1,dalidatpath_2
        character(len=100000) line
        character(len=6) key
        real zcut
c
        real score
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer idom1,idom2,ide,lali,nbest,ibest
        real zmax,rmsd,z,x,x1
        integer nres1_old,nres2_old

        read(line(1:n),*,err=11,end=99) key, cd1cd2,nblock,
     $          (l1(i),r1(i),i=1,nblock),(l2(i),r2(i),i=1,nblock)
11      continue
        !write(*,*) '#dowork_dp:line ',line(1:n)
!        write(*,*) '#dp in: ', key, cd1cd2,nblock,
!     $          (l1(i),r1(i),i=1,nblock),(l2(i),r2(i),i=1,nblock)
        cd1=cd1cd2(1:5)
        cd2=cd1cd2(6:10)
c
!       write(*,*) 'innerloop ',cd1,cd2,nres1,nres2,score,zcut
        if(cd1.ne.oldcd1) then
          oldcd1=cd1
          call dpsetup(cd1,ndom1,domns1,domseglist1,nres1,ca1,d1,
     $  90,dalidatpath_1,seq1)
!          write(*,*) '#dpsetup1 ',cd1,nres1,ndom1
c         file-not-found handle
          if(ndom1.eq.0) return
        end if
        if(cd2.ne.oldcd2) then
          oldcd2=cd2
          call dpsetup(cd2,ndom2,domns2,domseglist2,nres2,ca2,d2,
     $  90,dalidatpath_2,seq2)
!          write(*,*) '#dpsetup2 ',cd2,nres2,ndom2,nres1
c         file-not-found handle
          if(ndom2.eq.0) return
        end if

c       idom1 * idom2
c
        zmax=0.0
        x1=0.0
        rmsd=9.9
        ide=0
        nres1_old=nres1
        nres2_old=nres2
        do idom1=1,ndom1
                do idom2=1,ndom2
!                       write(*,*) 'dopair ',idom1,idom2
                        call dopair(idom1,idom2,
     $  l1,r1,l2,r2,nblock,cd1,cd2,z,x)
                        if(z.gt.zmax) zmax=z
                        if(idom1.eq.1.and.idom2.eq.1) x1=x
                end do
        end do
c
c       print protein pair if any domain pair had z>zcut
c
c        write(*,*) '#zmax',zmax,(zmax.ge.zcut),x1,' ',cd1,' ',cd2
        if(zmax.ge.zcut) then
!        write(*,*) '#l1dp',(l1(i),r1(i),i=1,nblock)
!        write(*,*) '#l2dp',(l2(i),r2(i),i=1,nblock)
                call compressblocks(nblock,l1,r1,l2,r2)
!                call getide(cd1,cd2,nblock,
!     $                  l1,r1,l2,r2,lali,ide,seq1,seq2,nres1)
               ! rmsd is filled by fssp program
!                call getrmsd(rmsd,nblock,l1,r1,l2,r2)
                write(outputunit,601) 'DCCP   1 ',x1,rmsd,abs(lali),
     $                  zmax,ide,nblock,cd1,cd2
                write(outputunit,605) 'alignment'
                write(outputunit,610) (l1(i),r1(i),i=1,nblock)
                write(outputunit,610) (l2(i),r2(i),i=1,nblock)
        end if

600     format(1x,a9,f8.1,f4.1,i4,16x,i4,3x,i4,16x,a5,1x,a5)
601     format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,3x,i4,16x,a5,1x,a5)
605     format(1x,a9)
610     format(8(i4,2x,i4))
620     format(a5)

!       write(*,*) 'innerloop done'

99      return
        end subroutine dowork_dp
c
c------------------------------------------------------------------------------
c
        subroutine getrmsd(rmsd,nblock,l1,r1,l2,r2)
        implicit none
        include 'parsizes.for'
        integer nblock,l1(nblock),r1(nblock),l2(nblock),r2(nblock)
        integer i,j,k,n,m,ierr,lali1,lali2
        real rmsd,x(3,maxres),y(3,maxres),w(maxres),ssq,u(3,3),t(3)

        n=0
        m=0
        lali1=0
        lali2=0
        do i=1,nblock
                write(*,*) '#getrmsd',i,l1(i),r1(i),l2(i),r2(i)
                if(r2(i).lt.l2(i)) r2(i)=l2(i)+r1(i)-l1(i)
                lali1=lali1+r1(i)-l1(i)+1
                lali2=lali2+r2(i)-l2(i)+1
                do j=l1(i),r1(i)
                        n=n+1
                        w(n)=1.0
                        do k=1,3
                                x(k,n)=ca1(k,j)
                        end do
                end do
                do j=l2(i),r2(i)
                        m=m+1
                        do k=1,3
                                x(k,m)=ca2(k,j)
                        end do
                end do
        end do
        if(n.ne.m) then
                write(*,*) '# getrmsd error:',n,m,lali1,lali2
                rmsd=-9.9
                return
        end if
        call u3b(w,x,y,n,0,ssq,u,t,ierr) ! calclate rms only
        write(*,*) '#u3b',n,ssq,ierr
        rmsd=sqrt(ssq/max(1,n))

        end subroutine getrmsd
c
c------------------------------------------------------------------------------
c
        subroutine getide(cd1,cd2,nblock,l1,r1,l2,r2,lali,intide,
     $  seq1,seq2,nres1)
        implicit none
        include 'parsizes.for'
        character*5 cd1,cd2
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres),lali
        integer intide
c
        character seq1(maxres),seq2(maxres),aligned(maxres)
        integer nres1,ali(maxres),i,j,k
        real xide
c
c        write(*,*) 'getide: ',cd1,cd2!,nblock,lali,intide,nres1

        do i=1,maxres
                ali(i)=0
        end do
        do i=1,nblock
                do j=l1(i),r1(i)
                        k=l2(i)+j-l1(i)
                        ali(j)=abs(k)
                end do
        end do
        ! aligned sequence string
        do i=1,maxres
                k=ali(i)
                if(k.eq.0) then
                        aligned(i)='.'
                else
                        aligned(i)=seq2(k)
                end if
        end do
        ! sequence identity
        xide=0.0
        lali=0
        do i=1,min(maxres,nres1)
                if(ali(i).ne.0) then
                  if(aligned(i).ge.'a'.and.seq1(i).ge.'a') then
                                xide=xide+100.0
                  else if(aligned(i).ge.'a'.and.seq1(i).eq.'C') then
                                xide=xide+100.0
                  else if(aligned(i).eq.'C'.and.seq1(i).ge.'a') then
                                xide=xide+100.0
                  else if(aligned(i).eq.seq1(i)) then
                                xide=xide+100.0
                  end if
                  lali=lali+1
                end if
        end do

        intide=xide/max(lali,1)

        return
        end subroutine getide
c
c------------------------------------------------------------------------------
c
        subroutine dopair(idom1,idom2,
     $    l1,r1,l2,r2,nblock,
     $    cd1,cd2,z,x)
c
c       extract idom1, idom2 from raw alignment
c       trim alignment
c       write out zscore,confined alignment
c
        implicit none
        include 'parsizes.for'
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer idom1,idom2
        character*5 cd1,cd2
c
        integer i,j,ali1(maxres),lali,len1,len2
        logical lactive1(maxres),lactive2(maxres)
        real z,x
c
c       mark domains
c
!       write(*,*) 'dopair',idom1,idom2,nres1,nres2,cd1,cd2
        do i=1,min(maxres,nres1)
                lactive1(i)=.false.
                ali1(i)=0
        end do
        x=0.0
        z=0.0
        len1=0
        do i=1,domns1(idom1)
                do j=domseglist1(i,1,idom1),domseglist1(i,2,idom1)
                        lactive1(j)=.true.
                        len1=len1+1
                end do
        end do
        do i=1,min(maxres,nres2)
                lactive2(i)=.false.
        end do
        len2=0
        do i=1,domns2(idom2)
                do j=domseglist2(i,1,idom2),domseglist2(i,2,idom2)
                        lactive2(j)=.true.
                        len2=len2+1
                end do
        end do
c
c       construct ali1
c
        lali=0
        do i=1,nblock
                do j=l1(i),r1(i)
                  if(lactive1(j).and.lactive2(l2(i)+j-l1(i))) then
                        ali1(j)=l2(i)+j-l1(i)
                        lali=lali+1
                  end if
                end do
        end do
!       write(*,*) 'lali',lali
c
        if(lali.eq.0) then
!               write(*,*) idom1,idom2,lali,' lali'
                return
        end if
c
c       calculate score
c
        x=totscore(ali1,nres1)
c        write(*,*) 'totscore =',x
        z=zscore(len1,len2,x)
c       write(*,*) '#z',z,x,' ',cd1,' ',cd2,nres1,nres2,len1,len2,lali
!
!       if(z.gt.0.90) then      !! write out all domain-domain pairs
                              !!      for calibration !
c               ! find blocks !
c        write(*,630) x,0.0,lali,l1(1),r1(nblock),l2(1),r2(nblock),
c     $                  -9,0,0,nblock,idom1,idom2,nres1,nres2,cd1,cd2
!               write(12,640) 1
!               write(12,650) (l1(i),r1(i),i=1,nblock)
!               write(12,650) (l2(i),r2(i),i=1,nblock)
!
!        end if
c
500     format('ZZZZ',f10.2,f10.1,3i5,2(2x,a5,i3))
630     format(' DCCP   1 ',f8.1,f4.1,6i4,2i2,3i3,2i4,2x,a5,1x,a5)
640     format(' alignment (resnos)',i10,':')
650     format(8(i4,' -',i4))
c
        return
        end subroutine dopair
c
c------------------------------------------------------------------------------
c
        function totscore(ali1,nres1) result(tots)
!        function totscore(ali1,d1,d2,nres1) result(tots)
        implicit none
        include 'parsizes.for'
        real tots
        integer nres1
        !integer*2 d1(maxres,maxres),d2(maxres,maxres)
        integer ali1(maxres)
        integer i,j,k,l,q,r,a(maxres),n
        real x
c
        tots=0.0
        n=0
        do i=1,min(maxres,nres1)
                if(ali1(i).ne.0) then
                        n=n+1
                        a(n)=i
                end if
        end do
        do i=1,n
                k=a(i)
                do j=1,n
                        l=a(j)
                        q=abs(ali1(k))
                        r=abs(ali1(l))
                        x=dpscorefun(d1(k,l),d2(q,r))
                        tots=tots+x
                end do
        end do
c
500     format(8(i4,' -',i4))
c
        return
        end function totscore
c
c----------------------------------------------------------------------
c
        function dpscorefun(a,b) result(s)
        implicit none
        include 'parsizes.for'
        real s
        integer*2 a,b
c
        real x,y,d0
        logical lela
        parameter(lela=.true.)
        parameter(d0=0.20)
c !!!   elastic uses weights !!!
        x=float(abs(a-b))/10
        if(lela) then
                y=float(a+b)/20
                if(y.gt.100) then
                        s=0.0
                else
                        if(y.gt.0) then
                          s=wght(nint(y))*(d0-x/y)
                        else
                          s=wght(nint(y))*d0
                        end if
                end if
        end if

        return
        end function dpscorefun
c
c-------------------------------------------------------------------------------
c
        subroutine dpweights
        implicit none
        include 'parsizes.for'
        integer i
c
        real x,enveloperadius
c
        enveloperadius=20.0
        x=1/(enveloperadius*enveloperadius)
        do i=0,100
                wght(i)=exp(-x*i*i)
        end do

        return
        end subroutine dpweights
c
c----------------------------------------------------------------------
c
        function distanceint2(a1,a2,a3,b1,b2,b3) result(d)
        implicit none
        real a1,a2,a3,b1,b2,b3
        integer*2 d

        d=nint(10.0*sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+
     $          (a3-b3)*(a3-b3)))

        return
        end function distanceint2
c
c-----------------------------------------------------------------------
c
        subroutine dpgetdist(ca,d,nres1)
        implicit none
        include 'parsizes.for'
        real ca(3,maxres)
        integer*2 d(maxres,maxres),x
        integer nres1
c
        integer i,j
c
        do i=1,min(maxres,nres1)
                d1(i,i)=0
                do j=1,i-1
                        x=distanceint2(ca(1,i),ca(2,i),ca(3,i),
     $                          ca(1,j),ca(2,j),ca(3,j))
                        d(i,j)=x
                        d(j,i)=x
                end do
        end do
c
        return
        end subroutine dpgetdist
c
c------------------------------------------------------------------------------
c
        function zscore(l1,l2,score) result(z)
        implicit none
        real z,score
        integer l1,l2
c
        real n12,mean,sigma,x
c
        n12=sqrt(float(l1*l2))
        x=min(n12,400.0)
        mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
        if(n12.gt.400.0) mean=mean+(n12-400.0)*1.0              ! hack !
        sigma=0.50*mean
        z=(score-mean)/max(1.0,sigma)

        return
        end function zscore
c
c------------------------------------------------------------------------------
c
        subroutine dpsetup(cd1,ndom,domns,domseglist,nres,ca,d,iunit,
     $          dalidatpath,seq)
        implicit none
        include 'parsizes.for'
        integer ndom,domns(maxdom),domseglist(maxseg,2,maxdom),iunit
        integer nres
        real ca(3,maxres)
        character*5 cd1
        integer*2 d(maxres,maxres)
        character seq(maxres)
c
        character*80 filnam,line,constructfilnam,dalidatpath
        character node_type(maxdom)
        integer i,j,k,idom,nseg
c
        nres=0
        ndom=0
!        write(*,*) 'setup ',cd1
        filnam=constructfilnam(cd1,dalidatpath,'.dat')
        open(iunit,file=filnam,status='old',err=19)
                i=0
200             read(iunit,530) line
!                write(*,*) '#read line ',line(1:len_trim(line))
                if(line(1:4).ne.'>>>>') goto 200
                read(line,550) ndom
                i=i+1
                if(i.eq.1) then
                        read(line(11:20),*) nres,nseg
                        do j=1,nseg
                                read(iunit,530) line
                        end do
                        read(iunit,*) ((ca(k,j),k=1,3),j=1,nres)
                end if
                if(i.lt.3) goto 200
                do idom=1,ndom
                        read(iunit,710,end=219) j,node_type(j),domns(j),
     $                          ((domseglist(i,k,j),k=1,2),i=1,domns(j))
                end do
c               -sequence
                read(iunit,540,end=219) (seq(i),i=1,nres)
219     close(iunit)
c
c       calculate distance matrix
c
        call dpgetdist(ca,d,nres)
c
c       keep only (+,*)-units
c
        ndom=0
        do i=1,j
          if(i.eq.1.or.node_type(i).eq.'+'.or.node_type(i).eq.'*') then
                        ndom=ndom+1
                        domns(ndom)=domns(i)
                        do k=1,domns(i)
                                domseglist(k,1,ndom)=domseglist(k,1,i)
                                domseglist(k,2,ndom)=domseglist(k,2,i)
                        end do
          end if
        end do
c
c       normal exit
c
        return
c
c       error exit
c
19      write(*,*) 'ERROR in fetchdomseglist: cannot open file',filnam
        ndom=0
c
530     format(a80)
540     format(10x,1x,64000a1)
550     format(10x,i5)
710     format(i4,1x,a1,13x,400i4)

        return
        end subroutine dpsetup
c
c------------------------------------------------------------------------------
c


        end module dp
c
c----------------------------------------------------------------------
c

