c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c	---------------------------------------------------------------------
c
	subroutine preparoutput(ali1i,ali1,blocksize,nres1,d1,d2,ca1,
     $		ca2,avscore,rms,alilen,tsil,tsil2)
c
c	calculates rms(popsize) and avscore(maxres,popsize) compressed
c	returns uncompressed ali1,avscore
c
	implicit none
	include 'gagasizes.for'
	integer nres1,blocksize,tsil(maxres),alilen(popsize)
        integer tsil2(maxres)
	real ca1(3,maxres),ca2(3,maxres)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer*2 ali1i(maxres,popsize),ali1(maxres0,popsize)
	real rms(popsize),avscore(maxres0,popsize),score(popsize)
c
c	local hack: summary reassigned in output()
c
	integer summary(11,popsize),i,a1(maxres)
c
	call uncompress(ali1i,ali1,avscore,tsil,nres1,tsil2)
	write(*,*) 'prepareoutput'
	call getrms(ali1,ca1,ca2,rms,blocksize,summary,nres1)
        write(*,*) '#rms, popsize',rms,popsize,lpretty
	do i=1,popsize
		alilen(i)=summary(1,i)
	end do
	if(lpretty) then
	  do i=1,popsize
		call getali0(ali1i,i,a1,blocksize,nres1)
		call getaveragescore(a1,i,d1,d2,nres1,avscore)
	  end do
	end if

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine uncompress(ali1i,ali1,avscore,tsil,nres1,tsil2)
	implicit none
	include 'gagasizes.for'
	integer*2 ali1(maxres0,popsize),ali1i(maxres,popsize)
	integer tsil(maxres),nres1,tsil2(maxres)
	real avscore(maxres0,popsize)
c
	integer i,ip,tmp(maxres0)
	real tmpa(maxres0)
c
	write(*,*) 'uncompress'
	do ip=1,popsize
		do i=1,maxres0
			tmp(i)=0
			tmpa(i)=0.0
		end do
		do i=1,nres1
		  if(tsil(i).gt.0.and.ali1i(i,ip).gt.0) then
			tmp(tsil(i))=tsil2(ali1i(i,ip))
			tmpa(tsil(i))=avscore(i,ip)
		  end if
		end do
		do i=1,maxres0
			ali1(i,ip)=tmp(i)
			avscore(i,ip)=tmpa(i)
		end do
	end do

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine output(code1,chainid1,code2,chainid2,nres1i,nres2i,
     $	seq1,seq2,struc1,struc2,resno1,resno2,blocksize,hdr1,hdr2,ca1,
     $	ca2,relacc1,relacc2,d1,d2,ali1i,score,nres0,tsil,nres02,tsil2,
     $  outunit)
	implicit none
	include 'gagasizes.for'
	character*4 code1,code2
	character chainid1,chainid2
	integer nres1,nres2,blocksize,relacc1(maxres0),relacc2(maxres0)
	character struc1(maxres0),struc2(maxres0),seq1(maxres0)
        character seq2(maxres0)
	character*4 resno1(maxres0),resno2(maxres0)
	character*80 hdr1,hdr2
	real ca1(3,maxres0),ca2(3,maxres0)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
        integer*2 ali1i(maxres,popsize)
	integer tsil(maxres),nres0,nres02,tsil2(maxres),nres1i,nres2i
        integer outunit
c
	character*80 pairfilnam
c
	integer i,j,k,l,ind,k1,alilen(popsize)
	integer d(popsize),ve(popsize)
	real score(popsize)
	integer*2 ali1(maxres0,popsize)
c
	integer summary(11,popsize)
	real rms(popsize)
	integer lcore1,lcore2,nide,lenali,l1,l2
	character aa,bb
c
	integer nblock(popsize),a1(maxres),b1(maxres,popsize),
     $		ib1(maxres,popsize),ib2(maxres,popsize),block(maxres)
	character id(60),se(60),st(60)
	character*80 filnam
	logical lali1
c
	real avscore(maxres0,popsize)
c
c	uncompress
c
	lverb=.false.
        ldebug=.false.
        lpretty=.false.
c        write(*,*) '# This is output ',(score(i),i=1,popsize)
c	write(*,*) code1,chainid1,code2,chainid2,nres1i,nres2i
c	write(*,*) nres0,nres02
c	write(*,300) (ali1i(i,1),i=1,nres1i)
	call preparoutput(ali1i,ali1,blocksize,nres1i,d1,d2,ca1,ca2,
     $		avscore,rms,alilen,tsil,tsil2)
	nres1=nres0
	nres2=nres02
c
c	summary table: 1 - nres, 2 - ifirst, 3 - ilast, 4 - jfirst, 5 - jlast
c		6 - %ide, 7 - nantiparallel, 8 - nreversal
c		9 - nblocks, 10 - %core1, 11 - %core2
c
	if(lverb)write(*,300) (ali1(i,1),i=1,nres1)
300	format(20i4)
	if(ldebug)write(*,*) ' dccp 1'
	do i=1,popsize
		summary(1,i)=alilen(i)
		do j=2,11
			summary(j,i)=0
		end do
	end do
	filnam=pairfilnam(code1,chainid1,code2,chainid2,'dccp')
	if(lverb)write(*,500) filnam
c !!!	print all output in one long file !!!
c	open(outunit,file=filnam,status='new')
c	write(outunit,510) code1,chainid1,hdr1
c	write(outunit,510) code2,chainid2,hdr2
c	call readali(pairfilnam(code1,chainid1,code2,chainid2,'ali '),
c     $		ali1,score)
	lcore1=0
!	do i=1,nres1
!		if(relacc1(i).le.corecutoff) lcore1=lcore1+1
!	end do
	lcore1=max(1,lcore1)
	lcore2=0
!	do i=1,nres2
!		if(relacc2(i).le.corecutoff) lcore2=lcore2+1
!	end do
	lcore2=max(1,lcore2)
	if(ldebug)write(*,*) ' dccp 2',(ali1(i,1),i=1,nres1)
c	call getsummary(ali1,summary)
	do i=1,popsize
		call getblocks(ali1,i,blocksize,block,k,nres1)
		summary(9,i)=k
		call getali0(ali1,i,a1,blocksize,nres1)
c		antiparallel segments
		k=0
		do j=1,nres1
			if((ali1(j,i).lt.0).and.(block(j).ne.k)) then
				k=block(j)
				summary(7,i)=summary(7,i)+1
			end if
		end do
		k=0
		l1=0
		l2=0
		nide=0
		lenali=0
		do j=1,nres1
		  if(a1(j).ne.0)then
			lenali=lenali+1
			if(a1(j).lt.k) summary(8,i)=summary(8,i)+1
			k=a1(j)
			if(summary(2,i).eq.0) summary(2,i)=j
			if(summary(3,i).lt.j) summary(3,i)=j
			k1=abs(a1(j))
			if(summary(4,i).eq.0) then
				summary(4,i)=k1
			else
				if(summary(4,i).gt.k1) summary(4,i)=k1
			end if
			if(summary(5,i).lt.k1) summary(5,i)=k1
!			if(relacc1(j).le.corecutoff) l1=l1+1
!			if(relacc2(k1).le.corecutoff) l2=l2+1
			aa=seq1(j)
			bb=seq2(k1)
			if((aa.eq.bb).or.((aa.ge.'a').and.(bb.ge.'a')))
     $				nide=nide+1
		  end if
		end do
		summary(10,i)=nint(100.0*float(l1)/float(lcore1))
		summary(11,i)=nint(100.0*float(l2)/float(lcore2))
		summary(6,i)=nint(100.0*float(nide)/float(max(1,lenali)))
		summary(1,i)=lenali
	if(ldebug)write(*,*) score(i),(summary(j,i),j=1,11),lenali,l1,
     $		lcore1,l2,lcore2
	end do
	if(ldebug)write(*,*) ' dccp 3'
c	call getrms(ali1,score,ca1,ca2,rms,blocksize,summary,nres1)
	if(ldebug)write(*,*) ' dccp 4'
c
c	set score to 0.0 if fewer than 4 blocks AND less than 30 residues
c
	call filter(score,summary,blocksize,4,min(nres1/2,nres2/2,30))
	if(ldebug)write(*,*) ' dccp 5'
	do i=1,popsize
		d(i)=i
		ve(i)=nint(score(i))
	end do
	if(ldebug)write(*,*) ' dccp 6'
	if(popsize.gt.1) call j_index_Qsort(1,popsize,popsize,d,ve)
c
c	exclude trivial submaxima 
c
	if(ldebug)write(*,*) ' dccp 7'
	do i=1,popsize
		call getali0(ali1,i,a1,blocksize,nres1)
		do j=1,maxres
			b1(j,i)=a1(j)
		end do
	end do
	do ind=1,popsize
	  i=d(ind)
	  if(score(i).gt.0.0)then
	    do j=popsize,ind+1,-1
	     if(score(d(j)).gt.0.0)then
	      do k=1,nres1
		if(b1(k,i).gt.0) then
		  if((b1(k,d(j)).gt.0).and.(b1(k,d(j)).eq.b1(k,i))) then
			score(i)=0.0
			goto 19
		  end if
		end if
	      end do
	     end if
	    end do
19	  end if
	end do
	call printsummary(score,rms,summary,d)
c	if(ldebug)write(*,*) ' dccp 8'
c	do ind=1,popsize
c	  i=d(popsize-ind+1)
c	  if(score(i).gt.0.0) write(outunit,505) ind,score(i),rms(i),
c     $	(summary(j,i),j=1,11),nres1,nres2,code1,chainid1,code2,chainid2
c	end do
c	call printblocks
	if(ldebug)write(*,*) ' dccp 9'
	do ind=1,printpopsize
	  i=d(popsize-ind+1)
	  if(score(i).gt.0.0) then
		call getblocks(ali1,i,blocksize,block,k,nres1)
		nblock(ind)=k
		do j=1,nblock(ind)
			ib1(j,ind)=0
			ib2(j,ind)=0
		end do
		do j=1,nres1
			if(block(j).gt.0)then
			  if(ib1(block(j),ind).eq.0) ib1(block(j),ind)=j
			  ib2(block(j),ind)=j
			end if
		end do
	  end if
	end do
	if(ldebug)write(*,*) ' dccp 10',(i,d(i),score(d(i)),i=1,popsize)
	do ind=1,printpopsize
	  i=d(popsize-ind+1)
	  if(score(i).gt.0.0) then
		call getali0(ali1,i,a1,blocksize,nres1)
c	write(*,530) ind
c	write(*,540) (ib1(j,ind),ib2(j,ind),j=1,nblock(ind))
c	write(*,540) (a1(ib1(j,ind)),a1(ib2(j,ind)),j=1,nblock(ind))

        write(outunit,*) 'WOLFITZ ',code1,chainid1,code2,chainid2,
     $          nblock(ind),(ib1(j,ind),ib2(j,ind),j=1,nblock(ind)),
     $         (a1(ib1(j,ind)),a1(ib2(j,ind)),j=1,nblock(ind))
!        write(outunit,*) 'DP',ind,score(i),rms(i),
!     $  (summary(j,i),j=1,11),nres1,nres2,code1,chainid1,code2,chainid2,
!     $  ind,(ib1(j,ind),ib2(j,ind),j=1,nblock(ind)),
!     $   (a1(ib1(j,ind)),a1(ib2(j,ind)),j=1,nblock(ind))
	  end if
	end do
	lali1=.false.
	if(ldebug)write(*,*) ' dccp 11'
!	do ind=1,printpopsize
!	  i=d(popsize-ind+1)
!	  if(score(i).gt.0.0) then
!		call getali0(ali1,i,a1,blocksize,nres1)
!		write(outunit,550) ind
!		write(outunit,560) (resno1(ib1(j,ind)),resno1(ib2(j,ind)),
!     $			j=1,nblock(ind))
!		write(outunit,560) (resno2(abs(a1(ib1(j,ind)))),
!     $			resno2(abs(a1(ib2(j,ind)))),j=1,nblock(ind))
!		lali1=.true.
!	   end if
!	 end do
c !!! no sequence output !
	if(lali.and.lali1)then
	  do i=1,nres1,60
		k=min(59,nres1-i)
		write(outunit,600) i,i+k,resno1(i),resno1(i+k)
		write(outunit,590) (struc1(i+j),j=0,k)
		write(outunit,570) code1,chainid1,(seq1(i+j),j=0,k)
		do ind=1,printpopsize
		  l=d(popsize-ind+1)
		  if(score(l).gt.0.0) then
			call getali0(ali1,l,a1,blocksize,nres1)
			do j=1,nres1
				a1(j)=abs(a1(j))
			end do
			do j=0,k
				se(j+1)=' '
				st(j+1)=' '
				id(j+1)=' '
				if(a1(i+j).ne.0) then
					se(j+1)=seq2(a1(i+j))
					st(j+1)=struc2(a1(i+j))
					if((struc1(i+j).ne.' ').and.
     $					  (st(j+1).eq.struc1(i+j)))
     $						id(j+1)=':'
					if((struc1(i+j).eq.'B').and.
     $					  (st(j+1).eq.'E')) id(j+1)=':'
					if((struc1(i+j).eq.'E').and.
     $					  (st(j+1).eq.'B')) id(j+1)=':'
					if(se(j+1).eq.seq1(i+j))
     $						id(j+1)='|'
					if((se(j+1).ge.'a').and.
     $						(seq1(i+j).ge.'a'))
     $						id(j+1)='|'
				end if
			end do
	write(outunit,590) (id(j),j=1,k+1)
	write(outunit,580) code2,chainid2,ind,(se(j),j=1,k+1)
	write(outunit,590) (st(j),j=1,k+1)
		  end if
		end do
	  end do
	end if
c
c	output pretty full alignment with boxes if best is without knots
c !!! no pretty !
	if(lpretty) then
	  i=d(popsize)
	  if((score(i).gt.0.0).and.(summary(8,i).eq.0)) then
		call pretty(a1,seq1,seq2,struc1,struc2,nres1,nres2,
     $			resno1,resno2,code1,code2,i,avscore,outunit)
	  end if
	end if
c	close(outunit)
500	format(/' generating output file ',a20)
505	format(' DCCP',i4,2x,f7.1,f4.1,6i4,2i2,3i3,2i4,1x,2(1x,a4,a1))
510	format(1x,a4,a1,4x,a70)
520	format(2i6,f10.3)
530	format(/' alignment (resnos)      ',i4,':')
540	format(8(i4,' -',i4))
550	format(/' alignment (PDB resnos)  ',i4,':')
560	format(8(a4,' -',a4))
570	format(1x,a4,a1,9x,60a1)
580	format(1x,a4,a1,i5,4x,60a1)
590	format(15x,60a1)
600	format(/' alignment residues ',i4,' -',i4,'  (',a4,' -',a4,')',
     $         /15x,6('         .'))

	return
	end
c
c	---------------------------------------------------------------------
c
        subroutine getali0(ali1,ind,a1,blocksize,nres1)
        implicit none
        include 'gagasizes.for'
        integer*2 ali1(maxres0,popsize)
        integer blocksize,a1(maxres),ind,nres1
c
        integer j,l,i1
 
        do j=1,min(maxres,nres1)
                a1(j)=0
        end do
        do j=1,min(maxres,nres1)
                i1=ali1(j,ind)
                if(i1.ne.0) then
                        do l=0,blocksize-1
                                a1(j+l)=i1+l
                        end do
                end if
        end do
 
        return
        end
c
c       ---------------------------------------------------------------------
c
	subroutine getaveragescore(a1,ind,d1,d2,nres1,avscore)
	implicit none
	include 'gagasizes.for'
	integer a1(maxres),ind,nres1
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real avscore(maxres0,popsize)
c
	integer i,j,n
	real s(maxres),smax(maxres),scorefun
	integer*2 a,b,ab
c
	do i=1,nres1
		s(i)=0.0
		smax(i)=0.0
	end do
	n=0
	do i=1,nres1
	 if(a1(i).ne.0) then
	  n=n+1
	  do j=1,nres1
	    if(a1(j).ne.0) then
		a=d1(i,j)
		b=d2(abs(a1(i)),abs(a1(j)))
		ab=(a+b)/2
		s(i)=s(i)+scorefun(a,b)
		smax(i)=smax(i)+scorefun(ab,ab)
	    end if
	  end do
	 end if
	end do
	do i=1,nres1
		avscore(i,ind)=s(i)/max(smax(i),1.0)
	end do

	return
	end
c
c	---------------------------------------------------------------------
c
	function scoresymbol(x)
	implicit none
	character scoresymbol
	real x
	character t(-1:9)
	data t/'-','0','1','2','3','4','5','6','7','8','9'/
	integer i

	i=int(10.0*x)
	if(i.lt.0) i=-1
	if(i.gt.9) i=9
	scoresymbol=t(i)

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine pretty(a1,seq1,seq2,struc1,struc2,nres1,nres2,
     $		resno1,resno2,code1,code2,ind,avscore,outunit)
	implicit none
	include 'gagasizes.for'
	integer a1(maxres),nres1,nres2,ind,outunit
	character struc1(maxres),struc2(maxres),seq1(maxres),seq2(maxres)
	character*4 resno1(maxres),resno2(maxres),code1,code2
	real avscore(maxres0,popsize)
c
	character line(2*maxres,5),scoresymbol
	character*4 r
	integer i0,i,j0,j,pos,p,k,x1(2*maxres),x2(2*maxres),i1,i2

	do i=1,2*maxres
		do k=1,5
			line(i,k)='.'
		end do
		x1(i)=0
		x2(i)=0
	end do
	i0=0
	j0=0
	pos=0
	do i=1,nres1
		if(a1(i).gt.0) then
			do p=i0+1,i-1
				line(pos+p-i0,1)=struc1(p)
				line(pos+p-i0,2)=seq1(p)
				x1(pos+p-i0)=p
			end do
			j=a1(i)
			do p=j0+1,j-1
				line(pos+p-j0,4)=seq2(p)
				line(pos+p-j0,5)=struc2(p)
				x2(pos+p-j0)=p
			end do
			do p=1,min(i-i0,j-j0)-1
				line(pos+p,3)='.'
			end do
			pos=pos+max(i-i0,j-j0)
			line(pos,1)=struc1(i)
			line(pos,2)=seq1(i)
			x1(pos)=i
			line(pos,3)=scoresymbol(avscore(i,ind))
			line(pos,4)=seq2(j)
			line(pos,5)=struc2(j)
			x2(pos)=j
			i0=i
			j0=a1(i)
		end if
	end do
	i=nres1
	do p=i0+1,i
		line(pos+p-i0,1)=struc1(p)
		line(pos+p-i0,2)=seq1(p)
	end do
	j=nres2
	do p=j0+1,j
		line(pos+p-j0,4)=seq2(p)
		line(pos+p-j0,5)=struc2(p)
	end do
	do p=pos+1,min(pos+nres2-j0,pos+nres1-i0)
		line(p,3)='.'
	end do
	pos=pos+max(i-i0,j-j0)
	write(outunit,*) ' '
	do i=1,pos,60
		k=min(59,pos-i)
		j=0
		i1=0
		do while((x1(i+j).eq.0).and.(j.lt.k))
			j=j+1
		end do
		if(j.le.k) i1=x1(i+j)
		j=0
		i2=0
		do while((x2(i+j).eq.0).and.(j.lt.k))
			j=j+1
		end do
		if(j.le.k) i2=x2(i+j)
		write(outunit,500) (line(i+j,1),j=0,k)
		if(i1.gt.0) r=resno1(i1)
		if(i1.eq.0) r='    '
		write(outunit,510) code1,r,(line(i+j,2),j=0,k)
		write(outunit,500) (line(i+j,3),j=0,k)
		if(i1.gt.0) r=resno2(i2)
		if(i1.eq.0) r='    '
		write(outunit,510) code2,r,(line(i+j,4),j=0,k)
		write(outunit,500) (line(i+j,5),j=0,k)
		write(outunit,*) ' '
	end do
500	format(10x,60a1)
510	format(1x,a4,a4,1x,60a1)

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine filter(score,summary,blocksize,minblocks,minlen)
	implicit none
	include 'gagasizes.for'
	real score(popsize)
	integer summary(11,popsize),blocksize,minblocks,minlen
c
	integer ind

	do ind=1,popsize
	  if((summary(1,ind).lt.minlen).and.(summary(9,ind).lt.minblocks))
     $		score(ind)=0.0
	end do

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine printsummary(score,rms,summary,d)
	implicit none
	include 'gagasizes.for'
	real score(popsize),rms(popsize)
	integer summary(11,popsize),d(popsize)
c
	integer i,j,ind

	do ind=1,popsize
	  i=d(popsize-ind+1)
	  write(*,500) ind,score(i),rms(i),
     $		(summary(j,i),j=1,10)
	end do
500	format('#summary',i4,7x,f7.1,f4.1,6i4,2i2,2i3)

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine getrms(ali1,ca1,ca2,rms,blocksize,summary,nres1)
c
c	fills rms(1..popsize)
c
	implicit none
	include 'gagasizes.for'
	integer blocksize,summary(11,popsize),nres1
	integer*2 ali1(maxres0,popsize)
	real ca1(3,maxres),ca2(3,maxres),rms(popsize)
c
	integer ind,a1(maxres),n,ierr,i,j
	real t(3),u(3,3),x(3,maxres),y(3,maxres),w(maxres),ssq

	do ind=1,popsize
		call getali0(ali1,ind,a1,blocksize,nres1)
		n=0
		do i=1,nres1
		  if(a1(i).ne.0)then
			n=n+1
			w(n)=1.0
			do j=1,3
				x(j,n)=ca1(j,i)
				y(j,n)=ca2(j,abs(a1(i)))
			end do
		  end if
		end do
		ssq=0.0
		if(n.gt.0) call u3b(w,x,y,n,0,ssq,u,t,ierr)
		rms(ind)=sqrt(ssq/max(n,1))
		summary(1,ind)=n
	end do

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine readali(filnam,ali1,score)
	implicit none
	include 'gagasizes.for'
	integer*2 ali1(maxres0,popsize)
	real score(popsize)
	character*80 filnam
c
	integer i,j,k,nres1,nres2
c
	open(90,file=filnam,status='old')
	read(90,500) nres1,nres2,k
	do i=1,popsize
		do j=1,maxres
			ali1(j,i)=0
		end do
	end do
	read(90,500) ((ali1(i,j),i=1,nres1),j=1,k)
	read(90,510) (score(i),i=1,k)
	close(90)
500	format(10i8)
510	format(8f10.2)

	return
	end
c
c	---------------------------------------------------------------------
c 
