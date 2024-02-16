c This module/program is part of DaliLite (c) L. Holm 1999
c
	program fssp
c
c   usage: ... fsspfilter.perl | pipeforwarder 1ppt dbsize dalidatpath
c		| fssp.$CPUARC | html.perl > 1ppt.html
c
c
	implicit none
	integer maxres
	parameter(maxres=6000)
c
c	step 1: read DCCP file -- assumed filtered, sorted
c	step 2: fetch DSSP info
c	step 3: statistics + topology
c	step 4: merge split blocks
c	step 5: write summary-line
c	step 6: write alignment-strings
c	step 7: write range-block
c	step 8: write u()t()-block
c
	character*5 cd1,cd2
	character*80 dalidatpath_1,dalidatpath_2,compnd1,compnd
	integer nprot2,nblock,l1(maxres),l2(maxres)
        integer r1(maxres),r2(maxres)
	real score,zscore
	integer i,j,k,ierr
c
	character seq(maxres),stru(maxres),aligned(maxres),seq1(maxres)
	character stru1(maxres),c
	integer nres,nchain,nres1,ali(maxres)
	real ca1(3,maxres),ca(3,maxres)
	character*5 resno1(0:maxres),resno(0:maxres)
	character*3 onetothree
	integer from1,from2,to1,to2
	character*6 strid1,strid2
	character*24 flags
	logical lpermut(maxres),lrevers(maxres),lincons(maxres)
	integer nrevers,npermut,lenali,ier,n,dbsize,ierr1
	real x(3,maxres),y(3,maxres),u(3,3),t(3),rmsd,ssq,w(maxres),ide
        real z(3,maxres)
	logical used(maxres)
	character*80 line,dssptemplate1,dssptemplate2,dssptemplate3

	character rotateaxis(3)
	real angle(3),x_deg,y_deg,z_deg
c
c	input parameters -- reads dccp file from stdin
c
        if(iargc().lt.4) stop "USAGE: fssp cd1 dbsize DAT1 DAT2"
        call getarg(1,cd1)
        call getarg(2,line)
        read(line,*) dbsize
        call getarg(3,dalidatpath_1)
        call getarg(4,dalidatpath_2)
c
c	step 0: write out dbsize
c
	write(*,470) dbsize
c
c	step 1: process query -- nchain, nres, cmpnd, seqlength
c
	call getdalidata(dalidatpath_1,cd1,99,nres1,seq1,stru1,compnd1,
     $		ca1,resno1)
	if(cd1(5:5).ne.' ') then
		strid1='code-X'
		strid1(6:6)=cd1(5:5)
	else
		strid1='code  '
	end if
	strid1(1:4)=cd1(1:4)
	write(*,480) strid1
	write(*,485) nres1,(resno1(i),i=1,nres1)
c
c	step 2: process DCCP file
c
	ierr=0
	nprot2=0
	do while (ierr.eq.0)
         call read_nextDCCP(cd1,cd2,score,zscore,nblock,l1,r1,l2,r2,
     $		ierr,5,maxres)
	 write(*,*) 'read_nextDCCP returned ',cd1,cd2,zscore,nblock,ierr
         if(ierr.eq.0) then
	  if(cd2(5:5).ne.' ') then
		strid2='code-X'
		strid2(6:6)=cd2(5:5)
	  else
		strid2='code  '
	  end if
	  strid2(1:4)=cd2(1:4)
          call getdalidata(dalidatpath_2,cd2,99,nres,seq,stru,compnd,
     $			ca,resno)
		write(*,*) '//'
		write(*,490) cd2
		write(*,521) '-compnd ','"',compnd,'"'
		do i=1,nres1
			ali(i)=0
		end do
		do i=1,nblock
			do j=l1(i),r1(i)
				k=l2(i)+j-l1(i)
				ali(j)=abs(k)
			end do
		end do
		! aligned sequence string, to compute pide 
		do i=1,nres1
			k=ali(i)
			if(k.eq.0) then
				aligned(i)='.'
			else
				aligned(i)=seq(k)
			end if
		end do
        write(*,530) '-sequen ','"',(seq(i),i=1,nres),'"' ! unaligned sequence
		! sequence identity
		ide=0.0
		do i=1,nres1
		  if(ali(i).ne.0) then
	if(aligned(i).ge.'a'.and.seq1(i).ge.'a') then
		ide=ide+100.0
	else if(aligned(i).ge.'a'.and.seq1(i).eq.'C') then
		ide=ide+100.0
	else if(aligned(i).eq.'C'.and.seq1(i).ge.'a') then
		ide=ide+100.0
	else if(aligned(i).eq.seq1(i)) then
		ide=ide+100.0
	end if
		  end if
		end do
        write(*,530) '-struct ','"',(stru(i),i=1,nres),'"' ! unaligned dssp-sequence
		! merge blocks
		do i=2,nblock
	if(l1(i).eq.r1(i-1)+1.and.l2(i).eq.r2(i-1)+1) then
		write(*,*) 'merge blocks'
		l1(i)=l1(i-1)
		l2(i)=l2(i-1)
		l1(i-1)=0
		l2(i-1)=0
	end if
		end do
		k=nblock
		nblock=0
		do i=1,k
			if(l1(i).ne.0) then
				nblock=nblock+1
				l1(nblock)=l1(i)
				r1(nblock)=r1(i)
				l2(nblock)=l2(i)
				r2(nblock)=r2(i)
			end if
		end do
		call new_topology(l1,r1,l2,r2,lrevers,lpermut,lincons,
     $			nrevers,npermut,nblock,lenali,maxres,used)
		! range lines
		do i=1,nblock
			from1=l1(i)
			to1=r1(i)
			from2=l2(i)
			to2=r2(i)
			flags="                        "
			if(lrevers(i)) flags(1:8)=  'REVERSED'
			if(lpermut(i)) flags(10:17)='PERMUTED'
			if(lincons(i)) flags(19:24)='DOUBLE'
            write(*,550) strid1,strid2,
     $				from1,to1,from2,to2,
     $				onetothree(seq1(from1)),resno1(from1),
     $				onetothree(seq1(to1)),resno1(to1),
     $				onetothree(seq(from2)),resno(from2),
     $				onetothree(seq(to2)),resno(to2),flags
		end do
		! summary line:
		! z,rmsd,lali,lseq2,%ide,revers,permut,nfrag,topo,compnd
		n=0
		do i=1,nres1
		  k=abs(ali(i))
		  if(k.ne.0) then
			n=n+1
			do j=1,3
				x(j,n)=ca(j,k) 	! superimpose second on query
				y(j,n)=ca1(j,i)
			end do
			w(n)=1.0
			!write(*,*) 'pair: ',n,i,k,(x(j,n),j=1,3),(y(j,n),j=1,3)
		  end if
		end do
		do i=1,3
			t(i)=0.0
			do j=1,3
				u(i,j)=0.0
			end do
		end do
		if(n.gt.0) then
			call u3b(w,x,y,n,1,ssq,u,t,ier)
			rmsd=sqrt(ssq/n)
			!write(*,*) 'u3b gives ',n,ssq,rmsd,ier,u,t
		else
			rmsd=0.0
		end if

		do i=1,3
            write(*,570)  strid1,strid2,i,(u(i,j),j=1,3),t(i)
		end do

                ! convert u() to rotate x|y|z
                call u2xyzrotate(u,rotateaxis,x_deg,y_deg,z_deg)
                if(rotateaxis(1).eq.'x') then
                        angle(1)=x_deg
                        angle(2)=y_deg
                        angle(3)=z_deg
                else
                        angle(1)=y_deg
                        angle(2)=z_deg
                        angle(3)=x_deg
                end if
                write(*,580) (rotateaxis(i),angle(i),i=1,3),
     $                  (t(i),i=1,3)
580     format('-transrotate',3(1x,a1,f10.4),1x,3f10.3)

		if(nrevers.eq.0.and.npermut.eq.0) then
			c='S'
		else
			c='N'
		end if
        write(*,560) strid1,strid2,zscore,rmsd,lenali,nres,
     $		  nint(ide/max(1,n)),nrevers,npermut,nblock,c,
     $		  compnd
	 end if
       end do

470	format('-dbsize ',i5)
480	format('-query ',a6)
485	format('-nres ',i5,64000(1x,a5))
486	format('-nchain ',i5)
490	format('-protein ',a5)
500	format(a5)
510	format(a80)
521	format(a8,1x,a1,a70,a1)
520	format(a8,1x,a1,a80,a1)
530	format(a8,1x,64000a1)
540	format(a8,1x,a1,a33,a1)
550	format('-ranges  "',a6,1x,a6,2x,i4,' -',i4,' <=> ',i4,' -',i4,
     $   '   (',a3,1x,a5,' - ',a3,1x,a5,' <=> ',a3,1x,a5,' - ',
     $   a3,1x,a5,')',2x,a24,'"')
560	format('-summar  "',a6,1x,a6,f5.1,f5.1,i5,i6,i5,i7,i7,i6,1x,a1,
     $   4x,a70,'"')
570	format('-matrix  "',a6,1x,a6,'  U(',i1,',.) ',3f10.6,f20.6,'"')

	return
	end
c
c-------------------------------------------------------------------------------
c
        subroutine read_nextDCCP(cd1,cd2,score,zscore,
     $          nblock,l1,r1,l2,r2,ierr,iunit,maxres)
        implicit none
c
c       reads dali-scores & Z-scores from domainparser output !
c
        character*5 cd1,cd2,c1,c2
        character*9 line
        real score,zscore
	integer maxres,nblock,l1(maxres),l2(maxres)
        integer r1(maxres),r2(maxres)
        integer i,ierr,iunit
        integer ide,lali
        real rmsd
	character*132 longline
c
cDCCP   1   1651.8 0.0 175    24.6         100      1                1hgeB 1hgeBc
100	read(iunit,500,end=999) longline
	read(longline,601,err=100) line,score,rmsd,lali,zscore,
     $          	ide,nblock,c1,c2
        if(line(1:9).eq.'DCCP   1 ') then
                if(c1.eq.cd1) then      ! expected order
10                      read(iunit,605,err=100,end=999) line
                        if(line(1:9).ne.'alignment') goto 10
                        read(iunit,610,err=100) (l1(i),r1(i),i=1,nblock)
                        read(iunit,610,err=100) (l2(i),r2(i),i=1,nblock)
                        cd2=c2
                else if(c2.eq.cd1) then ! reverse pair !
20                      read(iunit,605,err=100,end=999) line
                        if(line(1:9).ne.'alignment') goto 20
                        read(iunit,610,err=100) (l2(i),r2(i),i=1,nblock)
                        read(iunit,610,err=100) (l1(i),r1(i),i=1,nblock)
                        cd2=c1
                else
                        goto 100
                end if
        else
                goto 100
        end if

500	format(a132)
600     format(1x,a9,f8.1,f4.1,i4,16x,i4,4x,i3,6x,2i4,2x,a5,1x,a5)
601     format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,4x,i3,16x,a5,1x,a5)
605     format(1x,a9)
610     format(8(i4,2x,i4))

!       normal exit
        ierr=0
        return

!       error exit
999     ierr=-1
        return
	end
c
c-----------------------------------------------------------------------------
c
	subroutine new_topology(l1,r1,l2,r2,lrevers,lpermut,lincons,
     $		nrevers,npermut,nblock,lenali,maxres,used)
	implicit none
	integer maxres,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	integer nblock,lenali
	logical lpermut(maxres),lrevers(maxres),lincons(maxres),
     $          used(maxres)
	integer nrevers,npermut
c
	integer j,k,l,m

		lenali=0
		nrevers=0
		npermut=0
		do k=1,maxres
			used(k)=.false.
		end do
		do j=1,nblock
			lrevers(j)=.false.
			lpermut(j)=.false.
			lincons(j)=.false.
			if(l2(j).le.r2(j)) then
				l=l2(j)
				m=r2(j)
			else
				l=r2(j)
				m=l2(j)
			end if
			do k=l,m
				if(.not.used(k)) then
					used(k)=.true.
				else
					lincons(j)=.true.
				end if
			end do
			lenali=lenali+r1(j)-l1(j)+1
			if(l2(j).gt.r2(j)) then
				nrevers=nrevers+1
				lrevers(j)=.true.
			end if
			if(j.gt.1) then
				if(.not.lrevers(j-1).and.lrevers(j))
     $					lpermut(j)=.true.
				if((l2(j).lt.l2(j-1)).and.
     $	(.not.lrevers(j-1))) then
					npermut=npermut+1
					lpermut(j)=.true.
				end if
				if((l2(j).lt.l2(j-1)).and.
     $	(lrevers(j-1)).and.(.not.lrevers(j))) then
					npermut=npermut+1
					lpermut(j)=.true.
				end if
				if((l2(j).gt.l2(j-1)).and.
     $	(lrevers(j-1))) then
					npermut=npermut+1
					lpermut(j)=.true.
				end if
			end if
		end do

	return
	end

