c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
	subroutine getgagaprotein2(code,chainid,nres,d,blocksize,
     $		hash,lsearch)
	implicit none
	include 'gagasizes.for'
	character chainid
	character*4 code
	integer nres,blocksize
	integer*2 d(maxres,maxres)
c		
	integer nseg,segtype(maxres),segfirst(maxres),seglast(maxres)
	character*80 hashfilnam
	integer*2 hash(maxres,maxres)
	logical lsearch(maxres,maxres)
c
	integer i,j,m,n
	integer*2 psquare(maxres,maxres)
c
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5,psquare)
c
	if(lverb)write(*,500) hashfilnam(code,chainid)
	call hashblock(nres,hash,psquare,blocksize)
	call segments(psquare,nres,d,blocksize,nseg,segtype,
     $          segfirst,seglast)
	if(lverb)write(*,510) nseg,nres,code,chainid
c
	do i=1,nres
		do j=1,nres
			lsearch(i,j)=(psquare(i,j).le.p600)
		end do
	end do
c
	call partitionchain(psquare,segfirst,seglast,lsearch,blocksize,
     $		nseg,nres)
	if(lverb)then
		n=0
		m=0
		do i=1,nres-blocksize+1
			do j=i+blocksize,nres-blocksize+1
				if(lsearch(i,j)) m=m+1
			end do
		end do
		write(*,*) n,' T lhash',m,' lsearch'
	end if
c
500	format(/' generating lookup tables ',a20)
510	format(i10,' segments ',i10,' residues in ',a4,2x,a1)
c
	return 
	end
c
c----------------------------------------------------------------------
c
	subroutine getgagaprotein1(code,chainid,nres,d,blocksize,
     $		dist,distoverlay,hash)
	implicit none
	include 'gagasizes.for'
	character chainid
	character*4 code
	integer nres,blocksize
	integer*2 d(maxres,maxres)
	integer*2 hash(maxres,maxres),
     $		dist(2,0:maxd),distoverlay(2,maxres,maxres) 
c		
	character*80 hashfilnam
c
	integer i,j,m,n
	logical lhash(maxres,maxres)
	integer*2 psquare(maxres,maxres)
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1(1),lhash)
	equivalence (jnk1(maxres*maxres+1),psquare)
c
	if(lverb)write(*,500) hashfilnam(code,chainid)
	call gethash(d,nres,hash,blocksize)
	call hashblock(nres,hash,psquare,blocksize)
	if(lverb)write(*,510) nres,code,chainid
c
	do i=1,nres
		do j=1,nres
			lhash(i,j)=(psquare(i,j).le.p720)
			distoverlay(1,i,j)=0
			distoverlay(2,i,j)=0
		end do
	end do
	do i=0,maxd
		dist(1,i)=0
		dist(2,i)=0
	end do
c
	if(lverb)then
		n=0
		m=0
		do i=1,nres-blocksize+1
			do j=i+blocksize,nres-blocksize+1
				if(lhash(i,j)) n=n+1
			end do
		end do
		write(*,*) n,' T lhash'
	end if
c
	call hashhist(blocksize,psquare,nres,dist,distoverlay,lhash)
c
500	format(/' generating lookup tables ',a20)
510	format(i10,' residues in ',a4,2x,a1)
c
	return 
	end
c
c----------------------------------------------------------------------
c
	subroutine partitionchain(psquare,segfirst,seglast,lsearch,
     $		blocksize,nseg,nres)
	implicit none
	include 'gagasizes.for'
	integer*2 psquare(maxres,maxres)
	integer segfirst(maxres), seglast(maxres),blocksize
	integer nres,nseg
	logical lsearch(maxres,maxres)
c
	integer seglist(maxres),a(maxres),b(maxres),i,j,k,l,m,n
c
c	keep only one square per residue per segment (smallest psquare)
c
	m=0
	l=0
	do i=1,nseg
	  if(seglast(i)-segfirst(i).gt.1) then
		m=m+1
		seglist(m)=i
		l=l+seglast(i)-segfirst(i)+1
	  end if
	end do
	if(lverb)write(*,*) m,' long segments ',l,' residues'
	do i=1,m
		do j=segfirst(seglist(i)),seglast(seglist(i))
			do k=segfirst(seglist(i)),seglast(seglist(i))
				lsearch(j,k)=.false.
			end do
		end do
		do j=i+1,m
			do k=segfirst(seglist(i)),seglast(seglist(i))
				a(k)=9999
			end do
			do k=segfirst(seglist(j)),seglast(seglist(j))
				b(k)=9999
			end do
			do k=segfirst(seglist(i)),seglast(seglist(i))
				do l=segfirst(seglist(j)),seglast(seglist(j))
					if(psquare(k,l).lt.a(k))a(k)=k
					if(psquare(k,l).lt.b(k))b(k)=l
					lsearch(k,l)=.false.
				end do
			end do
			do k=segfirst(seglist(i)),seglast(seglist(i))
				lsearch(k,a(k))=.true.
			end do
			do k=segfirst(seglist(j)),seglast(seglist(j))
				lsearch(b(k),k)=.true.
			end do
		end do
	end do
	m=0
	n=0
	do i=1,nres-blocksize+1
		do j=i+1,nres-blocksize+1
			if(lsearch(i,j)) m=m+1
			n=n+1
		end do
	end do
	if(lverb)write(*,*) m,' squares to check of ',n,
     $		 ' save ',float(n-m)/n*100,' % '

	return 
	end
c
c----------------------------------------------------------------------
c
	subroutine segments(psquare,nres,d,blocksize,
     $		nseg,segtype,segfirst,seglast)
	implicit none
	include 'gagasizes.for'
	integer*2 d(maxres,maxres)
	integer nres,blocksize
	integer nseg,segtype(maxres),segfirst(maxres),seglast(maxres)
	integer*2 psquare(maxres,maxres)
c
	integer i,j,map(maxres),di(maxres),ve(maxres)
	logical lsearch(maxres,maxres),fok,l1
	real x,score50,scorefun
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1,lsearch)
c	
	score50=0.5*blocksize*blocksize*scorefun(0,0)
	do i=1,nres
	 	do j=1,nres
			lsearch(i,j)=.false.
		end do
		map(i)=i
	end do
c
c	quick segment definition as repeated backbone patterns
c
	do i=1,nres-blocksize+1
		do j=i+1,nres-blocksize+1
		  if((.not.lsearch(i,i)).and.(.not.lsearch(j,j)).and.
     $		     (abs(psquare(i,i)-psquare(j,j)).le.20)) then
			l1=fok(d,d,i,i,j,j,blocksize,x,0)
			if(l1.and.(x.gt.score50)) then
				lsearch(i,j)=.true.
				lsearch(j,i)=.true.
			end if
		  end if
		end do
	end do
c
c	cluster starting from most-neighbours fragment
c
	do i=1,nres
		map(i)=i
		ve(i)=0
		di(i)=i
	end do
	do i=1,nres
		do j=1,nres
			if(lsearch(i,j)) ve(i)=ve(i)+1
		end do
	end do
	call j_index_Qsort(1,nres,maxres,di,ve)
	do i=nres,1,-1
		do j=1,nres
			if(.not.(lsearch(j,j)).and.(lsearch(di(i),j))) then
				map(j)=di(i)
				lsearch(j,j)=.true.
			end if
		end do
	end do
	if(lverb) write(*,510) (map(i),i=1,nres)
510	format(20i4)
c	segments are square-blocks
	call getsegments(map,nres-blocksize+1,nseg,segtype,
     $          segfirst,seglast)

	return 
	end
c
c-----------------------------------------------------------------------------
c
	subroutine getsegments(defi,nres,nseg,segtype,segfirst,seglast)
	implicit none
	include 'gagasizes.for'
	integer nres,defi(nres)
	integer nseg,segtype(maxres),segfirst(maxres),seglast(maxres)
c
	integer i,j
c
	i=1
	j=1
	nseg=0
10	do while((i.lt.nres).and.(defi(i).eq.defi(j)))
		i=i+1
	end do
	if(defi(j).gt.0)then
		nseg=nseg+1
		segtype(nseg)=defi(j)
		segfirst(nseg)=j
		seglast(nseg)=i-1
	end if
	j=i
	if(i.lt.nres)goto 10
	if(lverb) write(*,500) nseg
	if(lverb)then
	  do i=1,nseg
		write(*,510) i,segtype(i),segfirst(i),seglast(i),
     $			seglast(i)-segfirst(i)+1
	  end do
	end if

500    format(i5,' segments')
510    format(i4,' type: ',i3,'  first: ',i5,'   last:  ',i5,
     $  '   length: ',i5)

	return
	end
c
c-----------------------------------------------------------------------------
c
	subroutine gethash(d1, nres1, hash1, blocksize)
	implicit none
	include 'gagasizes.for'
	integer nres1, blocksize
	integer*2 hash1(maxres,maxres)
	integer*2 d1(maxres,maxres)
	real x
	integer i,j,k
c
c	hash1 returns 1*blocksize sum-of-distances
c
	do i=1,nres1
		do j=1,nres1
			hash1(i,j)=0
		end do
	end do
	do i=1,nres1
		x=0.0
		do k=0,blocksize-1
			x=x+d1(i,1+k)/10
		end do
		hash1(i,1)=min(nint(x),maxd)
		do j=2,nres1-blocksize+1
			x=x-d1(i,j-1)/10+d1(i,j+blocksize-1)/10
			hash1(i,j)=min(nint(x),maxd)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------------=
c
	subroutine hashblock(nres,hash,psquare,blocksize)
	implicit none
	include 'gagasizes.for'
	integer nres,blocksize
	integer*2 hash(maxres,maxres), psquare(maxres,maxres)
c
	integer i,j,k,s
c
c	pblock returns the square-sum-of-distances 
c
	do i=1,nres
		do j=1,nres
			psquare(i,j)=0
		end do
	end do
	do i=1,nres-blocksize+1
		do j=i,nres-blocksize+1
			s=0
			do k=0,blocksize-1
				s=s+hash(i+k,j)
			end do
			psquare(i,j)=s
			psquare(j,i)=s
		end do
	end do
	
	return
	end
c
c----------------------------------------------------------------------------=
c
	subroutine hashhist(blocksize,psquare,nres,dist,distoverlay,lhash)
	implicit none
	include 'gagasizes.for'
	integer blocksize,nres
	integer*2 dist(2,0:maxd),distoverlay(2,maxres,maxres)
	integer*2 psquare(maxres,maxres)
	logical lhash(maxres,maxres)
c
	integer i,j,k,n,i1,j1
c
c	dist,distoverlay return lookups to square-sum-of-distances
c	                 initialized in getprotein
c
	n=0
	do i=1,nres-blocksize+1
		do j=i+blocksize,nres-blocksize+1
		  if(lhash(i,j)) then
			n=n+1
			k=min(psquare(i,j),maxd)
			if(dist(1,k).eq.0) then
				dist(1,k)=i
				dist(2,k)=j
			else
				i1=dist(1,k)
				j1=dist(2,k)
				distoverlay(1,i,j)=distoverlay(1,i1,j1)
				distoverlay(2,i,j)=distoverlay(2,i1,j1)
				distoverlay(1,i1,j1)=i
				distoverlay(2,i1,j1)=j
			end if
		  end if
		end do
	end do
	if(lverb)write(*,*) ' hashist: ',n
 
	return
	end
c
c----------------------------------------------------------------------------
c

