c This module/program is part of DaliLite (c) L. Holm 1999
c
	subroutine triplet(ngene,genepool,trip,n3)
	implicit none
	include 'gagasizes.for'
	integer ngene
	integer*2 genepool(5,ngene)
	integer trip(3,10000)
c
c	make list of all triplets in genepool A-B & A-C & B-C = A-B-C
c
	integer hash(maxres,-maxres:maxres),overlay(maxpair)
	integer nx,i1,i2,j1,j2,k,list(maxpair),i,j,a1,a2,ilist(maxpair)
	integer ii,ij,n3,n,ni,igene,jgene,kgene,x1,x2,b1,b2,c1,c2
	logical seedlist(maxpair)
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1,hash)
	integer jnk6(dim6)
	common /junk6/jnk6
	equivalence (jnk6(1),overlay)
	equivalence (jnk6(maxpair+1),list)
	equivalence (jnk6(maxpair+maxpair+1),ilist)
	equivalence (jnk6(maxpair+maxpair+maxpair+1),seedlist)
c
	do igene=1,ngene
		overlay(igene)=0
		seedlist(igene)=.false.
	end do
	do i=1,maxres
		do j=-maxres,maxres
			hash(i,j)=0
		end do
	end do
	nx=0
	n=0
	do igene=1,ngene
	  if(genepool(5,igene).gt.0) then
		i1=genepool(1,igene)
		i2=genepool(3,igene)
		j1=genepool(2,igene)
		j2=genepool(4,igene)
		if(i1.ge.j1) write(*,*) 'BUG',igene,i1,i2,j1,j2
		k=hash(i1,i2)
		if(k.eq.0) then
			hash(i1,i2)=igene
			n=n+1
			list(n)=igene
		else
			overlay(igene)=overlay(k)
			overlay(k)=igene
			nx=nx+1
		end if
	  end if
	end do
	if(lverb)write(*,*) ngene,' genes in hash table,',
     $ nx,' in overlay,',n,' occupied'
c
c	combinatorial search
c
	n3=0
	do i=1,n
		igene=list(i)
		ni=0
		do while(igene.gt.0)
			ni=ni+1
			ilist(ni)=igene
			igene=overlay(igene)
		end do
		do ii=1,ni
			igene=ilist(ii)
			a1=genepool(1,igene)
			a2=genepool(3,igene)
			b1=genepool(2,igene)
			b2=genepool(4,igene)
			do ij=ii+1,ni
			  jgene=ilist(ij)
			  c1=genepool(2,jgene)
			  c2=genepool(4,jgene)
c			  case b1<c1
			  if(b1.lt.c1) then
				kgene=hash(b1,b2)
				do while(kgene.gt.0)
					x1=genepool(2,kgene)
					x2=genepool(4,kgene)
					if(x1.ne.c1) goto 19
					if(x2.ne.c2) goto 19
	if(ldebug) write(*,*) 'triplet 1',n3,a1,a2,b1,b2,c1,c2,
     $		igene,jgene,kgene
	seedlist(igene)=.true.
	seedlist(jgene)=.true.
	seedlist(kgene)=.true.
					n3=n3+1
					trip(1,n3)=igene
					trip(2,n3)=jgene
					trip(3,n3)=kgene
19					kgene=overlay(kgene)
				end do
c			  case b1>=c1
			  else if(b1.ge.c1) then
				kgene=hash(c1,c2)
				do while(kgene.gt.0)
					x1=genepool(2,kgene)
					x2=genepool(4,kgene)
					if(x1.ne.b1) goto 29
					if(x2.ne.b2) goto 29
	if(ldebug) write(*,*) 'triplet 2',n3,a1,a2,b1,b2,c1,c2,
     $		igene,jgene,kgene
	seedlist(igene)=.true.
	seedlist(jgene)=.true.
	seedlist(kgene)=.true.
					n3=n3+1
					trip(1,n3)=igene
					trip(2,n3)=jgene
					trip(3,n3)=kgene
29					kgene=overlay(kgene)
				end do
			  end if
			end do
		end do
	end do
	n=0
	do i=1,ngene
		if(seedlist(i)) n=n+1
	end do
	if(lverb)write(*,*) n,' genes in seedlist'

	return
	end
		
