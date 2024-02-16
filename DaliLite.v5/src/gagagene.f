c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
	subroutine getgenes(code1,chainid1,code2,chainid2,d1,d2,nres1,
     $	nres2,blocksize,ngene,genepool,lpreali,preali1)
	implicit none
	include 'gagasizes.for'
	character chainid1,chainid2
	character*4 code1,code2
	integer nres1,nres2,blocksize
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
	logical lpreali
c	 
	character*80 pairfilnam,hashfilnam
c
	character*80 filnam
	integer*2 dist1(2,0:maxd),distoverlay1(2,maxres,maxres),
     $		hash1(maxres,maxres)
	integer*2 dist2(2,0:maxd),distoverlay2(2,maxres,maxres),
     $		hash2(maxres,maxres) 
	integer*2 genepool(5,maxpair)
	integer ngene,ngene2
c
	logical lhash1(maxres,maxres),lhash2(maxres,maxres),
     $		test1(maxres,-maxres:maxres),test2(maxres,-maxres:-1)
	logical fok,x,lself
	integer i,j,a3,a4,b
	real x0,scorefun
	integer*2 genepool2(5,maxpair),buffer(5,maxpair*2)
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1,test1)
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23(1),distoverlay1)
	equivalence (jnk23(maxres*maxres+1),distoverlay2)
	integer jnk4(dim4)
	common /junk4/jnk4
	equivalence (jnk4(1),hash1)
	equivalence (jnk4(maxres*maxres/2+1),hash2)
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5,test2)
	integer jnk6(dim6)
	common /junk6/jnk6
	equivalence (jnk6,lhash1)
	equivalence (jnk6,lhash2)
	equivalence (jnk6,buffer)

        integer*2 zero
        parameter(zero=0)
c
c	hash proteins
c
	call getgagaprotein1(code1,chainid1,nres1,d1,blocksize,
     $		dist1,distoverlay1,hash1)
	call getgagaprotein2(code1,chainid1,nres1,d1,blocksize,
     $		hash1,lhash1)
	call getgagaprotein1(code2,chainid2,nres2,d2,blocksize,
     $		dist2,distoverlay2,hash2)
c !!!	normalized cutoff
	filnam=pairfilnam(code1,chainid1,code2,chainid2,'gene')
	if(lverb)write(*,500) filnam
c
c !!! symmetry: partition A check all of B; partition B check rest of A
c
	if(.not.ltube) then
	  b=blocksize-1
	  do i=1,nres1
		do j=-nres2,nres2
			test1(i,j)=.false.
		end do
	  end do
	  do i=1,nres1-b
		do j=-nres2,-b-1
			test1(i,j)=fok(d1,d2,i,i,-j,-j,blocksize,x0,3)
		end do
		do j=1,nres2-b
			test1(i,j)=fok(d1,d2,i,i,j,j,blocksize,x0,0)
 		end do
	  end do
c	  do i=1,nres1
c		write(*,600) (test1(i,j),j=-nres2,nres2)
c	  end do
c
c	  suppress 0-diagonal in comparing a protein with itself
c
	  lself=.false.
	  if(hashfilnam(code1,chainid1).eq.hashfilnam(code2,chainid2))then
		lself=.true.
		write(*,*) ' identical proteins, suppressing 0...10-diagonals'
		do i=1,nres1
			do j=i,i+10
				test1(i,j)=.false.
				test1(j,i)=.false.
			end do
		end do
	  end if
	else if(lpreali) then
c
c	  special case: limit search space to tube around imported prealignment
c
		b=blocksize-1
		do i=1,nres1
			do j=-nres2,nres2
				test1(i,j)=.false.
			end do
		end do
		do i=1,nres1
			if(preali1(i).gt.0) test1(i,preali1(i))=.true.
		end do
c
c	  special case: ltube=.true. => use main diagonal +- 20 residues only
c
	else
		write(*,*) ' tube option: suppressing all but +-20 diagonal '
		do i=1,nres1
			do j=-nres2,nres2
				test1(i,j)=.false.
			end do
			do j=max(1,i-20),min(nres2,i+20)
				test1(i,j)=fok(d1,d2,i,i,j,j,blocksize,x0,0)
			end do
		end do
	end if
c
c	two passes, equal chances for symmetry !
c
c	first call A-B
	ngene=0
	call getpairs(blocksize,dist2,distoverlay2,
     $		test1,hash1,hash2,d1,d2,ngene,genepool,
     $		lhash1,dist1,distoverlay1,lpreali)
	if(lverb)write(*,*) ngene, ' genes a-B'
c	second call B-A
	if(lverb)write(*,*) 'test2'
	do i=1,max(nres1,nres2)
		do j=1,i
			x=test1(i,j)
			test1(i,j)=test1(j,i)
			test1(j,i)=x
		end do
	end do
	do i=1,nres2
		do j=-nres1,-blocksize
			test2(i,j)=test1(-j-b,-i-b)
		end do
	end do
	do i=1,nres2
		do j=-nres1,-blocksize
			test1(i,j)=test2(i,j)
		end do
	end do
c	do i=1,nres2
c		write(*,600) (test1(i,j),j=-nres1,nres1)
c	end do
	call getgagaprotein2(code2,chainid2,nres2,d2,blocksize,hash2,
     $          lhash2)
	ngene2=0
	if(.not.lself)call getpairs(blocksize,dist1,distoverlay1,
     $		test1,hash2,hash1,d2,d1,ngene2,genepool2,
     $		lhash2,dist2,distoverlay2,lpreali)
	if(lverb)write(*,*) ngene2, ' genes b-A'
c	transpose B-A genes to A-B !
	do i=1,ngene2
		a3=1
		a4=1
		if (genepool2(3,i).lt.0) then
			buffer(1,i)=-genepool2(3,i)-b
			buffer(3,i)=-genepool2(1,i)-b
		else
			buffer(1,i)=genepool2(3,i)
			buffer(3,i)=genepool2(1,i)
		end if
		if (genepool2(4,i).lt.0) then
			buffer(2,i)=-genepool2(4,i)-b
			buffer(4,i)=-genepool2(2,i)-b
		else
			buffer(2,i)=genepool2(4,i)
			buffer(4,i)=genepool2(2,i)
		end if
		buffer(5,i)=genepool2(5,i)
	end do
c	append genepool to buffer, compress in sortpurge
	do i=ngene2+1,ngene2+ngene
		do j=1,5
			buffer(j,i)=genepool(j,i-ngene2)
		end do
	end do
	ngene=ngene2+ngene
	call purgelist(buffer,ngene)
	call sortpurge(buffer,genepool,ngene,
     $		blocksize*blocksize*scorefun(zero,zero))	
c
c	
c

c
	if(lsavegenes) call savegenepool(filnam,ngene,genepool)
500	format(/' generating gene file ',a20)
600	format(80l1)
c
	return 
	end
c
c----------------------------------------------------------------------------=
c
	subroutine getpairs(blocksize,dist2,distoverlay2,test1,
     $		hash1,hash2,d1,d2,n1,w,lhash1,dist1,distoverlay1,lpreali)
	implicit none
	include 'gagasizes.for'
	integer blocksize
	integer*2 w(5,maxpair)
	integer*2 dist2(2,0:maxd),distoverlay2(2,maxres,maxres),
     $		hash1(maxres,maxres),hash2(maxres,maxres) 
	integer*2 dist1(2,0:maxd),distoverlay1(2,maxres,maxres)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer n1
	logical test1(maxres,-maxres:maxres),lhash1(maxres,maxres),lpreali
c
	integer maxpair2
	parameter(maxpair2=maxpair*2)
	integer*2 buffer(5,maxpair2)
	logical fok
	integer i,j,k,i2,j2,p,q,nx,b,ix,jx,ntest,mtest
	real x0,score50,scorefun,rowcut(0:6),colcut(0:6)
        integer*2 zero
        parameter(zero=0)

	score50=0.5*blocksize*blocksize*scorefun(zero,zero)
c
c	2-D search for superposable blocks in CA-CA distance maps
c
	nx=0
	b=blocksize-1
c
c	check closest contacts until maxpairs reached !
c
	ntest=0
	mtest=0
	do p=1,plimit
	  if(lverb)write(*,*)p
	  i=dist1(1,p)
	  j=dist1(2,p)
	  do while (i.gt.0) 
	    if(.not.lhash1(i,j)) goto 9
	    do k=0,blocksize-1
		rowcut(k)=0.20*hash1(i+k,j)
		colcut(k)=0.20*hash1(j+k,i)
	    end do
c
c 	    look around at 20 % radius of p !!!
c
	    do k=max(1,int(0.81*p)),min(maxd,int(1.223*p))
		i2=dist2(1,k)
		j2=dist2(2,k)
		do while (i2.gt.0) 
			ntest=ntest+1
c!!! prealignment
			if(lpreali) then
				if(.not.test1(i,i2)) goto 19
c				if(.not.test1(j,j2)) goto 19
				if((i.gt.j).and.(i2.le.j2)) goto 19
				if((i.lt.j).and.(i2.ge.j2)) goto 19
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+q,j2)).gt.rowcut(q)) goto 19
	 if(abs(hash1(j+q,i)-hash2(j2+q,i2)).gt.colcut(q)) goto 19
				end do
				if(fok(d1,d2,i,j,i2,j2,blocksize,x0,0))
     $					call putgene(i,j,i2,j2,0,x0,
     $	buffer,n1,blocksize)
c !!! topology switch
			else if(ltop) then
				if(.not.test1(i,i2)) goto 19
				if(.not.test1(j,j2)) goto 19
				if((i.gt.j).and.(i2.le.j2)) goto 19
				if((i.lt.j).and.(i2.ge.j2)) goto 19
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+q,j2)).gt.rowcut(q)) goto 19
	 if(abs(hash1(j+q,i)-hash2(j2+q,i2)).gt.colcut(q)) goto 19
				end do
				if(fok(d1,d2,i,j,i2,j2,blocksize,x0,0))
     $					call putgene(i,j,i2,j2,0,x0,
     $	buffer,n1,blocksize)
			else if(.not.lpara) then
c
c !!! ltop=.false. & lpara=.false. --> chain reversal disallowed
c
				if(.not.test1(i,i2)) goto 19
				if(.not.test1(j,j2)) goto 19
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+q,j2)).gt.rowcut(q)) goto 19
	 if(abs(hash1(j+q,i)-hash2(j2+q,i2)).gt.colcut(q)) goto 19
				end do 
				if(fok(d1,d2,i,j,i2,j2,blocksize,x0,0))
     $					call putgene(i,j,i2,j2,0,x0,
     $	buffer,n1,blocksize)
			
			else
c
c !!! ltop=.false. & lpara=.true. --> chain reversal allowed
c ++ 0 +- 1 -+ 2 -- 3
c
			if((.not.test1(i,i2)).and.(.not.test1(i,-i2))) goto 19
		 	if((.not.test1(j,j2)).and.(.not.test1(j,-j2))) goto 19
c	case 0 ++
				if(.not.test1(i,i2)) goto 119
				if(.not.test1(j,j2)) goto 119
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+q,j2)).gt.rowcut(q)) goto 119
	 if(abs(hash1(j+q,i)-hash2(j2+q,i2)).gt.colcut(q)) goto 119
				end do 
				if(fok(d1,d2,i,j,i2,j2,blocksize,x0,0))
     $					call putgene(i,j,i2,j2,0,x0,
     $	buffer,n1,blocksize)
c	case 1 +-
119				if(.not.test1(i,i2)) goto 219
				if(.not.test1(j,-j2)) goto 219
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+q,j2)).gt.rowcut(q)) goto 219
	 if(abs(hash1(j+q,i)-hash2(j2+b-q,i2)).gt.colcut(q)) goto 219
				end do 
				if(fok(d1,d2,i,j,i2,j2+b,blocksize,x0,1))
     $					call putgene(i,j,i2,j2,1,x0,
     $	buffer,n1,blocksize)
c	case 2 -+
219				if(.not.test1(i,-i2)) goto 319
				if(.not.test1(j,j2)) goto 319
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+b-q,j2)).gt.rowcut(q)) goto 319
	 if(abs(hash1(j+q,i)-hash2(j2+q,i2)).gt.colcut(q)) goto 319
				end do 
				if(fok(d1,d2,i,j,i2+b,j2,blocksize,x0,2))
     $					call putgene(i,j,i2,j2,2,x0,
     $	buffer,n1,blocksize)
c	case 3 --
319				if(.not.test1(i,-i2)) goto 19
				if(.not.test1(j,-j2)) goto 19
				do q=0,blocksize-1
	 if(abs(hash1(i+q,j)-hash2(i2+b-q,j2)).gt.rowcut(q)) goto 19
	 if(abs(hash1(j+q,i)-hash2(j2+b-q,i2)).gt.colcut(q)) goto 19
				end do 
				if(fok(d1,d2,i,j,i2+b,j2+b,blocksize,x0,3))
     $					call putgene(i,j,i2,j2,3,x0,
     $	buffer,n1,blocksize)
			end if
			if(n1.ge.maxpair2-100) goto 99
			mtest=mtest+1
19			ix=distoverlay2(1,i2,j2)
			jx=distoverlay2(2,i2,j2)
			i2=ix
			j2=jx
		end do
	    end do
9	    ix=distoverlay1(1,i,j)
	    jx=distoverlay1(2,i,j)
	    i=ix
	    j=jx
	  end do
	end do
99	if(lverb)write(*,*) ' end of gene generation at p=',p
	write(*,*) ' ntest,mtest ',ntest,mtest
c 	skip sort if n1<maxpair -- but must convert from buffer to w !
	call sortpurge(buffer,w,n1,2*score50)
	if(lverb)write(*,500) n1,nx

500	format(i10,' genes saved ',i10,' failed foks ')

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine putgene(ai,ak,ai2,aj2,xtype,x,buffer,n1,blocksize)
	implicit none
	include 'gagasizes.for'
	integer ai,ak,ai2,aj2,xtype,blocksize
	real x
	integer maxpair2,n1
	parameter(maxpair2=maxpair*2)
	integer*2 buffer(5,maxpair2)
c
	integer b,c
c
c !!! quit if buffer is full !!!
c
	  if(n1.eq.maxpair2) return
	if(x.lt.0.5) return
c
	  if((xtype.eq.0).or.(xtype.eq.1)) then
		b=ai2
	  else
		b=-(ai2+blocksize-1)
c 		error !
c		if(b.gt.-blocksize) return
	  end if
	  if((xtype.eq.0).or.(xtype.eq.2))then
		c=aj2
	  else
		c=-(aj2+blocksize-1)
c 		error !
c		if(c.gt.-blocksize) return
	  end if
c	  write(*,*) ' putgene: ',n1,ai,ak,b,c,xtype,x
	  n1=n1+1
	  buffer(1,n1)=ai
	  buffer(2,n1)=ak
	  buffer(3,n1)=b
	  buffer(4,n1)=c
	  buffer(5,n1)=nint(x)

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine sortpurge(buffer,w,n1,maxscore)
c
c	returns maxpair best genes in w and buffer and sets n1 to maxpair
c
	implicit none
	include 'gagasizes.for'
	integer maxpair2
	parameter(maxpair2=maxpair*2)
	integer n1
	integer*2 buffer(5,maxpair2)
	integer*2 w(5,maxpair)
	real maxscore
c
	integer i,n,p,q,ve(maxpair2),hist(0:100),nx
c
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5,ve)
c
c	return maxpair best genes from buffer in w
c
	if(lverb)write(*,*) n1, ' in to sortpurge'
	nx=0
	if(n1.le.maxpair) then
		do i=1,n1
			do p=1,5
				w(p,i)=buffer(p,i)
			end do
		end do
	else
c		divide scores in 100 bins, find maxpair-limit from histogram
		do q=0,100
			hist(q)=0
		end do
		do i=1,n1
			q=min(100,nint(100.0*buffer(5,i)/maxscore))
			if(q.lt.0) then
				q=0
				nx=nx+1
			end if
			ve(i)=q
			hist(q)=hist(q)+1
		end do
		n=0
		q=100
		do while((q.gt.0).and.(n.lt.maxpair))
			n=n+hist(q)
			q=q-1
		end do
		q=q+1
		if(lverb)write(*,*) ' %score cutoff: ',q, ' histogram:'
		if(lverb)write(*,*) (hist(i),i=0,100)
		n=0
		do i=1,n1
			if(ve(i).gt.q) then
				n=n+1
				do p=1,5
					w(p,n)=buffer(p,i)
				end do
			end if
		end do	
		n1=n
	end if
	if(lverb)write(*,*) n1,' out from sortpurge'

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine savegenepool(filnam,ngene,genepool)
	implicit none
	include 'gagasizes.for'
	integer ngene
	character*80 filnam
	integer*2 genepool(5,maxpair)
	integer i,j

c	call genehash(ngene,genepool,hashi,hashj,di,dj)
	open(90,file=filnam,status='new')
	write(90,500) ngene,((genepool(j,i),j=1,5),i=1,ngene)
	close(90)
500	format(10i8)

	return
	end
c
c----------------------------------------------------------------------------=
c
	function fok(d1,d2,i1,j1,i2,j2,blocksize,x,xtype)
	implicit none
	include 'gagasizes.for'
	logical fok
	integer*2 d1(maxres,maxres), d2(maxres,maxres)
	real x
	integer i1,j1,i2,j2,blocksize
	integer i,j,k1,k2,xtype
	real scorefun
c
c
c
	x=0.0
	if(xtype.eq.0) goto 10
	if(xtype.eq.1) goto 20
	if(xtype.eq.2) goto 30
	if(xtype.eq.3) goto 40
	  write(*,*) ' unknown xtype in fok ',xtype
	  fok=.false.
	  return
10	  do i=0,blocksize-1
		k1=i1+i
		k2=i2+i
		do j=0,blocksize-1
			x=x+scorefun(d1(k1,j1+j),d2(k2,j2+j))
			if(x.lt.fokcutoff) goto 99
		end do
	  end do
	goto 99
20 	 do i=0,blocksize-1
		k1=i1+i
		k2=i2+i
		do j=0,blocksize-1
			x=x+scorefun(d1(k1,j1+j),d2(k2,j2-j))
			if(x.lt.fokcutoff) goto 99
		end do
	  end do
	goto 99
30	  do i=0,blocksize-1
		k1=i1+i
		k2=i2-i
		do j=0,blocksize-1
			x=x+scorefun(d1(k1,j1+j),d2(k2,j2+j))
			if(x.lt.fokcutoff) goto 99
		end do
	  end do
	goto 99
40	  do i=0,blocksize-1
		k1=i1+i
		k2=i2-i
		do j=0,blocksize-1
			x=x+scorefun(d1(k1,j1+j),d2(k2,j2-j))
			if(x.lt.fokcutoff) goto 99
		end do
	  end do

99	fok=(x.gt.0.0) 

	return
	end
c
c----------------------------------------------------------------------------=
c
	subroutine purgelist(buffer,ngene)
c
c	delete identical genes from buffer by setting score to zero
c
	implicit none
	include 'gagasizes.for'
	integer maxpair2,ngene
	parameter(maxpair2=maxpair*2)
	integer*2 buffer(5,maxpair2)
c
	integer test(maxres,-maxres:maxres),tmp(maxres,-maxres:maxres)
	integer overlay(maxpair2)
	integer i1,j1,i2,j2,k1,k2,i,j,n,nx
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1,test)
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23,tmp)
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5,overlay)
c
c	initialize
c
	do i=1,maxres
		do j=-maxres,maxres
			test(i,j)=0
			tmp(i,j)=0
		end do
	end do
c
c	sort genelist (i1-j1,i2-j2) so that i1<j1
c
	do i=1,ngene
		overlay(i)=0
		i1=buffer(1,i)
		j1=buffer(2,i)
		i2=buffer(3,i)
		j2=buffer(4,i)
		if(i1.gt.j1) then
			buffer(1,i)=i2
			buffer(2,i)=j2
			buffer(3,i)=i1
			buffer(4,i)=j1
		end if
	end do
c
c	fill test and overlay
c
	do i=1,ngene
		i1=buffer(1,i)
		j1=buffer(2,i)
		i2=buffer(3,i)
		j2=buffer(4,i)
		j=test(i1,i2)
		if(j.eq.0) then
			test(i1,i2)=i
		else
			overlay(i)=overlay(j)
			overlay(j)=i
		end if
	end do
c
c	check for multiple copies
c
	n=0
	nx=0
	do i=1,ngene
		i1=buffer(1,i)
		j1=buffer(2,i)
		i2=buffer(3,i)
		j2=buffer(4,i)
		j=test(i1,i2)
		n=n+1
		tmp(j1,j2)=n
c
c		only start chains from the gene in test()
c
	  if(j.eq.i) then
		j=overlay(j)
		do while(j.gt.0) 
			k1=buffer(2,j)
			k2=buffer(4,j)
c
c			test() says (i1-i2) match, annul if (j1-j2) match 
c
			if(tmp(k1,k2).eq.n) then
				buffer(5,j)=0
				nx=nx+1
			else
				tmp(k1,k2)=n
			end if
			j=overlay(j)
		end do
	  end if
	end do
	if(lverb)write(*,500) nx,ngene

500	format(' scores of ',i5,' redundant copies set to zero',i10,
     $ ' genes ')

	return
	end

c
c----------------------------------------------------------------------------=
c

