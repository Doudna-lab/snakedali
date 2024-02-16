c This module/program is part of DaliLite (c) L. Holm 1999
c
	program selfdccp
	implicit none
	include 'parsizes.for'
        integer i,nprot2,nres
        character*5 list2(maxprot),cd1
        character*80 constructfilnam,listfile,filnam,dalidatpath
c
c	write(*,*) 'this program reads list2 and looks up nres in dalidatpath'
c	write(*,*) 'and writes out self-alignment in DCCP format'
c	write(*,*) 'continue with "dp dalidatpath DCCP 2.0"'
        if(iargc().lt.1) stop "USAGE: selfdccp DALIDATDIR < cdlist"
	call getarg(1,dalidatpath)
c
        call getlist(listfile,5,list2,nprot2)
	do i=1,nprot2
		cd1=list2(i)
                if(cd1(5:5).eq.' ') cd1(5:5)='_'
                filnam=constructfilnam(cd1,dalidatpath,'.dat')
                open(99,file=filnam,status='old',err=19)
		read(99,500) nres
		if(cd1(5:5).eq.'_') cd1(5:5)=' '
                write(*,511) cd1,cd1,1,1,nres,1,nres
!		write(*,510) nres,cd1,cd1,1,nres,1,nres
19		close(99)
	end do
c
500	format(10x,i5)
510	format(' DCCP   1    -99.9 0.0',i4,'    99.9         100',
     $  '      1                ',a5,1x,a5,/' alignment',2(/i4,2x,i4))
511     format('WOLFITZ ',2a5,5i6)
520	format(a80)
c
	end
!DCCP   1    210.7 1.9  51     6.8          20      4                1fxd  1fdx
!alignment
!  1    18  20    24  29    41  44    58
!  1    18  21    25  26    38  39    53
c
c----------------------------------------------------------------------
c
        subroutine getlist(filnam,iunit,list,nprot)
        implicit none
        include 'parsizes.for'
        character*5 list(maxprot)
        integer nprot,iunit
        character*(*) filnam
c
        character*5 cd
c
        nprot=0
        if(iunit.ne.5) open(iunit,file=filnam,status='old')
10      read(iunit,500,end=19) cd
        if(nprot.eq.maxprot) then
c               write(*,*) 'WARNING: skip reading list after maxprot',maxprot
                goto 19
        end if
        nprot=nprot+1
        list(nprot)=cd
        goto 10
19      if(iunit.ne.5) close(iunit)
c
500     format(a5)
c

        return
        end
c
c----------------------------------------------------------------------
c
