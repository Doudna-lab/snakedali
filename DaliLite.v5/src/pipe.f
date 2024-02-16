c This module/program is part of DaliLite (c) L. Holm 1999
c
        program refinput
        implicit none
        include 'parsizes.for'
        integer idom,i,j,k,ali(maxres)
        integer s,n,ndom
        integer domns(maxdom),domseglist(maxseg,2,maxdom)
        character*80 filnam,line,constructfilnam
        character*80 dalidatpath_1,dalidatpath_2
        character*5 cd1,cd2,oldcd1,cd,cdx
        integer a(2,maxres),b(2,maxres),nseg,keep(maxres)
        logical l,lfirst
        character cdom
        integer itox(maxdom),oldidom,ix
c
c
c
        if(iargc().lt.2) stop "USAGE: pipe dalidatpath_1 dalidatpath_2"
        call getarg(1,dalidatpath_1)
        call getarg(2,dalidatpath_2)
        oldcd1='     '
        oldidom=0
        lfirst=.true.
10      read(*,610,end=19,err=10) cd1,cd2,idom,s,nseg,
     $    (a(1,i),a(2,i),i=1,nseg),(b(1,i),b(2,i),i=1,nseg)
        l=(oldcd1.eq.cd1.and.oldidom.eq.idom)
        if(.not.l) then
                l=.true.
                oldcd1=cd1
                cd=cd1
                if(cd(5:5).eq.' ') cd(5:5)='_'
c
c               output
c
                if(.not.lfirst) write(*,500) 'END  '
                lfirst=.false.
                if(cd1(5:5).eq.'_') cd1(5:5)=' '
                cdom=' '
                !!!!if(idom.gt.1.or.idom.eq.1) cdom='*'
                write(*,500) cd1,' ',cdom
c
c               domain ranges for first protein
c
c               write(*,*) 'begin domain',idom,cdom
c                if(cdom.eq.'*') then
c                  write(*,*) domns(idom)
c        write(*,*) ((domseglist(i,k,idom),k=1,2),i=1,domns(idom))
c                end if
c               write(*,*) '...done'
        end if
        if(cd2(5:5).eq.'_') cd2(5:5)=' '
        write(*,500) cd2,'*'
        n=0
        do i=1,nseg
                if(a(1,i).gt.0) then
                        n=n+1
                        keep(n)=i
                end if
        end do
        write(*,*) n
        write(*,*) (a(1,keep(i)),a(2,keep(i)),i=1,n)
        write(*,*) (b(1,keep(i)),b(2,keep(i)),i=1,n)
        goto 10
19      close(98)
        write(*,500) 'END  '
c        write(*,500) 'END  '

500     format(a5,a1,a1)
510     format(6x,2a5,i4,i20,12i4,8(/8x,18i4))
520     format(15x,i5)
530     format(a80)
540     format(10x,i10)
551     format(10x,i5)
550     format(54x,100i4)
610     format(6x,2a5,i4,i20,i4,64000i4)
700     format(a5,2x,2i4)
710     format(i4,15x,400i4)

        end
c
c-------------------------------------------------------------------------------c
c
