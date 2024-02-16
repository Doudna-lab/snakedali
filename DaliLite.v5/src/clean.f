c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c-------------------------------------------------------------------------------
c	
	subroutine mclean(imax0,ngene,genepool,gblocksize,blocksize,
     $		nres1,nres2,d1,d2,bestscore,lsave,initfragali,bestali1,
     $		itrimx,ikillx,bestresali1,score,killcutx,ind,lkill,
     $		lfrozen,outfragali,reflen,nclonex,expectscore)
	implicit none
	include 'gagasizes.for'
	integer imax,ngene,blocksize,gblocksize,nres1,nres2,itrimx,reflen
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real totscore,expectscore
	logical lsave
	integer*2 genepool(5,maxpair),initfragali(maxres,popsize),
     $		bestali1(maxres),outfragali(maxres)
	integer nclonex
c
c	input
c		imax0 -- max # Monte Carlo steps
c		itrimx -- trim every itrimx'th step
c		gblocksize -- used by genepool
c		blocksize -- blocksize used by alignment
c			N.B.: chops up genepool to genes of length blocksize 
c		initfragali -- set of blocksize-genes to start from
c		ltop -- TRUE means sequential constraint
c		lsave -- TRUE means use rescore,dscore; FALSE means calculate all
c		ikillx -- call killer every ikillx'th step
c		bestresali1 -- saved alignments for killer
c		score -- scores of saved alignments
c		killcutx -- max allowed sequence identity for killer
c		ind -- identifies alignment for killer
c		lfrozen -- for killer
cc		reflen -- reference lenali for killer (in last refinement)
c		nclonex -- # clones to test in initfragali
c		expectscore -- return after 1st cycle if totscore is below
c	output
c		outfragali -- optimized set of genes
c		bestscore -- score of initfragali
c		bestali1 -- per residue alignment
c		lkill -- TRUE if alignment is converged (return if so)
c		lfrozen -- TRUE if alignment is converged
c
	integer cnt(maxres),prev1(0:maxres+1),next1(0:maxres+1)
	integer nextres(0:maxres+1),next2(0:maxres+1)
	integer*2 ali1(maxres),ali2(maxres),fragali1(maxres)
        integer*2 fragali2(maxres)
 	logical lset(maxpair),ltest(maxres,maxres)
	integer ncand, candij(2,maxpair*4)
	integer ntest,testij(2,maxpair),test(maxres,-maxres:maxres)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer ikillx,ind,imax0
	real score(popsize),killcutx
	integer*2 bestresali1(maxres,popsize)
	logical lkill,lfrozen(popsize)
c
c	local bookkeeping
c		cnt(i) -- multiplicity of genes at residue i
c		ali1(i),ali2(j) -- per residue alignment
c		fragali1(i),fragali2(j) -- per gene alignment
c		prev1(i),next1(i) -- lookup with respect to fragali1(i)
c		next2(j) -- lookup with respect to fragali2(j)
c		nextres(i) -- lookup with respect to ali1(i)
c		lset(ix) -- TRUE if gene is in fragali1
c		ltest(i,j) -- TRUE if residue pair has been appended to cand-list
c		ncand -- length of candidate (residue pairs i,j) list
c		candij(2,ix) -- pointer from candidate to residue pair
c		rescore(i,j) -- worth of residue pair in current alignment
c		dscore(i,j) -- how much residue pair would add to current ali
c		ntest -- length of active gene list
c		testij(2,it) -- pointer from gene list to first residue pair
c		test(i,j) -- pointer from genes first residue pair to gene list
c
c	make sure to delete incompatible residues before adding new ones !
c
c
c	functions
c
	real getp
	real ran
c
c	locals
c
	integer istep,nacc,ndel,ncand0,lenali
	logical lfr
	real alfa,delta,bestscore
	real dx,p
	integer i,j,k,ix
	integer nrem,rem(2,maxres),nnew,new(2,maxres),ngenedel
	integer di(maxpair),ve(maxpair),q,imp,iclone,genedel(maxres)
	real impscore
c
	integer jnk1(dim1)
	common /junk1/jnk1
	equivalence (jnk1,test)
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23(1),rescore)
	equivalence (jnk23(maxres*maxres+1),dscore)
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5(1),candij)
	equivalence (jnk5(8*maxpair+1),testij)
	equivalence (jnk5(10*maxpair+1),di)
	equivalence (jnk5(11*maxpair+1),ve)
	equivalence (jnk5(12*maxpair+1),lset)
	integer jnk6(dim6)
	common /junk6/jnk6
	equivalence (jnk6,ltest)
c
	imax=imax0
	if(lverb)write(*,*)'This is Monte Carlo',imax,ngene,blocksize,
     $		gblocksize,ltop,lsave
c
c	initialize
c
	bestscore=0.0
	call initialize(cnt,ali1,ali2,fragali1,fragali2,prev1,next1,
     $		ltest,ncand,rescore,dscore,nres1,nres2,ntest,nextres,
     $		test,next2)
c
c	load initali
c
	lenali=0
	totscore=0.0
	do i=1,nres1
		j=initfragali(i,1)
		if(j.ne.0) then
			call appendcand(i,j,blocksize,
     $				lset,ltest,dscore,ncand,
     $				candij,ali1,ali2,nres1,d1,d2,
     $				nextres,lsave,ntest,test,testij)
c			call printcandlist(ncand,candij,dscore,rescore)
			call testaddition(i,j,lsave,blocksize,prev1,
     $	next1,ali1,cnt,nres1,d1,d2,rescore,dscore,
     $	dx,nrem,rem,nnew,new,ngenedel,genedel,
     $	nres2,fragali1,fragali2,nextres,next2)

			call changeconfig(ngenedel,genedel,nrem,rem,
     $	i,j,nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)

		end if
	end do
	if(lverb)write(*,*) 'initial totscore',totscore
c!!!	call gettotscore(ali1,d1,d2,nres1,totscore)
c	write(*,500) (initfragali(i,1),i=1,nres1)
c	write(*,500) (ali1(i),i=1,nres1)
c
c	Monte Carlo step
c
       do iclone=1,nclonex
	lfr=.false.
	ncand0=0
	alfa=50.0
	delta=0.0
	imp=0
	istep=0
	impscore=0.0
	do while(istep.lt.imax)
		istep=istep+1
c
c		kill !
c
		if(mod(istep,ikillx).eq.1) then
			call killer(bestresali1,ali1,bestscore/100,
     $	nres1,blocksize,score,killcutx,lfrozen,lkill,
     $	ind,reflen)
			if(lkill) then
	if(lverb)write(*,*) ' killing converged ',ind
				return
			end if
		end if
c
c		heatshock if score improvement is stalled
c
		if(istep-imp.eq.1) then
			alfa=50.0
			delta=0.0
		else if(istep-imp.eq.2) then
			alfa=1.0/sqrt(max(10000.0,totscore))
			delta=alfa/1000
		end if
		ndel=0
		nacc=0
		alfa=alfa+delta
c
c		pull more genes as long as find more candidates
c
		if(.not.lfr) then
			call pullgene(gblocksize,blocksize,ngene,
     $	genepool,lset,ltest,dscore,ncand,candij,fragali1,next1,
     $	ali1,ali2,nres1,nres2,d1,d2,nextres,
     $	lsave,ntest,test,testij)
c			call printcandlist(ncand,candij,dscore,rescore)
			if(ncand.eq.ncand0) then
				lfr=.true.
			else
				imax=imax+1
			end if
			ncand0=ncand
		end if
c
c		adding phase  WARM
c
c
c		test in randomized order
c
		do ix=1,ntest
			di(ix)=ix
			ve(ix)=nint(100.0*ran(seed))
		end do
		call j_index_Qsort(1,ntest,maxpair,di,ve)
		do q=1,ntest
			ix=di(q)
			i=testij(1,ix)
			j=testij(2,ix)
			if(fragali1(i).ne.j) then
				if(fragali1(i).eq.j) then
	write(*,*) 'WARNING: duplicate gene',i,j
					goto 19
     				end if
				call testaddition(i,j,lsave,blocksize,
     $	prev1,next1,ali1,cnt,nres1,d1,d2,
     $	rescore,dscore,dx,nrem,rem,nnew,new,
     $	ngenedel,genedel,nres2,
     $	fragali1,fragali2,nextres,next2)
				p=getp(alfa*dx)
				if(ran(seed).lt.p) then
c				  write(*,*)'accepted',i,j,dx

				  call changeconfig(ngenedel,genedel,
     $		nrem,rem,i,j,
     $		nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $		nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $		cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $		nres1,outfragali,bestali1,dx)
				  if(totscore.gt.impscore) then
					imp=istep
					impscore=totscore
				  end if
				end if
19			end if
		end do
c!!!		call gettotscore(ali1,d1,d2,nres1,totscore)
c
c		trimming phase COLD
c
		if(mod(istep,itrimx).eq.0) then
		 i=0
		 do while(next1(i).le.nres1)
		  i=next1(i)
		  j=fragali1(i)
c		  Trim from ends only
		  if((cnt(i).eq.1).or.(cnt(i+blocksize-1).eq.1))then
		ngenedel=1
		genedel(1)=i
		call testdeletion(ngenedel,genedel,blocksize,cnt,
     $			rescore,ali1,nextres,nres1,d1,d2,lsave,dx,
     $			nrem,rem)
c			p=getp(alfa*dx)
c			if(ran(seed).lt.p) then
c			only accept improvements !
			if(dx.gt.0.0) then
c				write(*,*)'trimming',i,j,dx

				  call changeconfig(ngenedel,genedel,
     $	nrem,rem,0,0,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)
				  if(totscore.gt.impscore) then
					imp=istep
					impscore=totscore
				  end if

			end if
		  end if
		 end do
		end if
		if(lverb)write(*,510) istep,imp,nacc,ndel,lenali,ncand,
     $			totscore,bestscore,alfa
c!!!		call gettotscore(ali1,d1,d2,nres1,totscore)
c		if(lsave)call print(ncand,candij,dscore,rescore,ali1,
c     $			next1,d1,d2,nres1)
c
c		compare to expected score
c
		if(totscore.lt.expectscore) return
	end do
c
c	mutate fragali1 to initfragali(iclone+1)
c
	if(iclone.lt.nclonex) then
	 if(lverb)write(*,*) ' clone',iclone
	 if(lverb)write(*,500) (fragali1(i),i=1,nres1)
	 if(lverb)write(*,500) (initfragali(i,iclone+1),i=1,nres1)
c
c		additions, replacements
c
		do i=1,nres1
		  j=initfragali(i,iclone+1)
		  k=fragali1(i)
		  if((j.ne.0).and.(j.ne.k)) then
				call testaddition(i,j,lsave,blocksize,
     $	prev1,next1,ali1,cnt,nres1,d1,d2,
     $	rescore,dscore,dx,nrem,rem,nnew,new,
     $	ngenedel,genedel,nres2,
     $	fragali1,fragali2,nextres,next2)
c				write(*,*) 'add',i,j
c				write(*,500) (rem(1,ix),rem(2,ix),ix=1,nrem)
c				write(*,500) (new(1,ix),new(2,ix),ix=1,nnew)

				call changeconfig(ngenedel,genedel,
     $	nrem,rem,i,j,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)

		  end if
		end do
c!!!		call gettotscore(ali1,d1,d2,nres1,totscore)
c
c		deletions
c
		do i=1,nres1
		  j=initfragali(i,iclone+1)
		  k=fragali1(i)
		  if((k.ne.0).and.(k.ne.j)) then
			ngenedel=1
			genedel(1)=i
			call testdeletion(ngenedel,genedel,blocksize,
     $		cnt,rescore,ali1,nextres,nres1,d1,d2,lsave,dx,
     $		nrem,rem)
c			write(*,*) 'delete',i,k
c			write(*,500) (rem(1,ix),rem(2,ix),ix=1,nrem)

			call changeconfig(ngenedel,genedel,
     $		nrem,rem,0,0,
     $		nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $		nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $		cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $		nres1,outfragali,bestali1,dx)

		  end if
		end do
c!!!		call gettotscore(ali1,d1,d2,nres1,totscore)
	end if
	imax=imax0
       end do
c
500	format(20i4)
510	format('step',6i6,2f10.2,f10.5)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine changeconfig(ngenedel,genedel,nrem,rem,i,j,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,lset,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,test,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)
c
c	input:
c		ngenedel,genedel -- genes to delete
c		i,j -- gene to add
c		nrem,rem -- residues to delete
c		nnew,new -- residues to add
c		istep -- for imp
c		dx -- change in totscore
c
c	output:
c		imp -- updated if score improved
c		totscore -- updated
c		outfragali -- bestever fragali1 updated
c		bestali1 -- bestever ali1 updated
c
	implicit none
	include 'gagasizes.for'
	integer ngenedel,genedel(maxres),nrem,rem(2,maxres),i,j,nnew
	integer new(2,maxres)
	integer*2 fragali1(maxres),fragali2(maxres),ali1(maxres)
        integer*2 ali2(maxres)
	integer next1(0:maxres+1),prev1(0:maxres+1),next2(0:maxres+1)
	integer nextres(0:maxres+1),candij(2,maxpair*4),ncand,lenali
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real rescore(maxres,maxres)
	real dscore(maxres,maxres),totscore,bestscore,dx
	integer cnt(maxres),ndel,nacc,test(maxres,-maxres:maxres)
	logical lset(maxpair),lsave
	integer nres1,blocksize
	integer*2 outfragali(maxres),bestali1(maxres)
c
	integer k
c
	call dodeletions(ngenedel,genedel,nrem,rem,fragali1,fragali2,
     $	ali1,ali2,next1,prev1,lset,candij,ncand,d1,d2,rescore,dscore,
     $		lenali,cnt,ndel,lsave, test,nextres,blocksize,next2)
	if(i.ne.0) call doadditions(i,j,nnew,new,fragali1,fragali2,
     $	ali1,ali2,next1,prev1,lset,candij,ncand,d1,d2,rescore,dscore,
     $	lenali,cnt,blocksize,nacc,lsave,test,nextres,next2,nres1)
	totscore=totscore+dx
	if(ldebug) then
		write(*,*) 'totscore now',totscore,nacc,ndel
		write(*,500) (ali1(k),k=1,nres1)
		write(*,500) (fragali1(k),k=1,nres1)
		write(*,500) (cnt(k),k=1,nres1)
		write(*,500) (nextres(k),k=1,nres1)
c		if(lsave)call print(ncand,candij,dscore,
c     $			rescore,ali1,next1,d1,d2,nres1)
		write(*,500) (rem(1,k),rem(2,k),k=1,nrem)
		write(*,500) (new(1,k),new(2,k),k=1,nnew)
		call gettotscore(ali1,d1,d2,nres1,totscore)
		call checkcnt(cnt,fragali1,nres1,blocksize)
	end if
	if(totscore.gt.bestscore) then
		call saveali(nres1,fragali1,outfragali,totscore,bestscore)
		call saveali(nres1,ali1,bestali1,totscore,bestscore)
	end if
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine print(ncand,candij,dscore,rescore,
     $		ali1,next1,d1,d2,nres1)
	implicit none
	include 'gagasizes.for'
	integer ncand,candij(2,maxpair*4),next1(0:maxres+1),nres1
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real dscore(maxres,maxres)
	real rescore(maxres,maxres)
	integer*2 ali1(maxres)
c
	real getrescore
c
	integer k,ix,jx
	real x,a,b
c
	write(*,*) 'This is print',ncand
	do k=1,ncand
		ix=candij(1,k)
		jx=candij(2,k)
		x=getrescore(ix,jx,d1,d2,ali1,next1,nres1)
		a=dscore(ix,jx)
		b=rescore(ix,jx)
		if(abs(a+b-x).gt.0.1) write(*,*)'dscore',ix,jx,a,b,x,jx-ali1(ix)
	end do
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine dodeletions(ngenedel,genedel,nrem,rem,fragali1,
     $	fragali2,ali1,ali2,next1,prev1,lset,candij,ncand,d1,d2,
     $	rescore,dscore,lenali,cnt,ndel,lsave,test,nextres,blocksize,
     $	next2)
	implicit none
	include 'gagasizes.for'
	integer ngenedel,genedel(maxres),nrem,rem(2,maxres)
	integer next2(0:maxres+1)
	integer*2 fragali1(maxres),fragali2(maxres),ali1(maxres)
        integer*2 ali2(maxres)
	integer next1(0:maxres+1),prev1(0:maxres+1),cnt(maxres),blocksize
	logical lset(maxpair)
	integer candij(2,maxpair*4),ncand
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer lenali,ndel,test(maxres,-maxres:maxres)
        integer nextres(0:maxres+1)
	logical lsave
c
	real scorefun
c
	integer irem,i,j,j1,k,l,ix,jx,in,ip
	real s,dx
c
	if(ldebug)then
		write(*,*) 'This is dodeletions',ngenedel
		write(*,500) (genedel(i),i=1,ngenedel)
		write(*,500) (rem(1,i),rem(2,i),i=1,nrem)
		write(*,500) (ali1(i),i=1,38)
		write(*,500) (nextres(i),i=0,39)
	end if
c
c	delete genes
c
	do irem=1,ngenedel
		i=genedel(irem)
		j=fragali1(i)
		fragali1(i)=0
		fragali2(abs(j))=0
		lset(test(i,j))=.false.
		do l=0,blocksize-1
			cnt(i+l)=cnt(i+l)-1
			if(cnt(i+l).lt.0) write(*,*) 
     $				' WARNING: subzero cnt ' ,i,j,i+l,j+l
		end do
c
c		update prev1,next
c
		ip=prev1(i)
		in=next1(i)
		do ix=i,in
			prev1(ix)=ip
		end do
		do ix=ip,i
			next1(ix)=in
		end do
		j1=abs(j)
		in=next2(j1)
		ix=j1
		do while((ix.gt.0).and.(next2(ix-1).eq.j1))
			ix=ix-1
			next2(ix)=in
		end do
	end do
c
c	delete residues
c
	do irem=1,nrem
		i=rem(1,irem)
		j=rem(2,irem)
c
c		update nextres
c
		in=nextres(i)
		ix=i
c		write(*,*) 'nextres',in,i,nextres(ix-1)
		do while((ix.gt.0).and.(nextres(ix-1).eq.i))
			ix=ix-1
			nextres(ix)=in
c			if(ldebug)write(*,*) 'set nextres',ix,' to ',in
		end do
		if(lsave) then
c
c			update dscore(x,y), rescore(x,y)
c	
			do k=1,ncand
				ix=candij(1,k)
				jx=candij(2,k)
				dx=0.0
				if((i.ne.ix).and.(j.ne.jx)) then
					s=-scorefun(d1(ix,i),d2(jx,j))
					dx=s+s
				end if
				if(ali1(ix).eq.jx) then
					rescore(ix,jx)=rescore(ix,jx)+dx
				else
					dscore(ix,jx)=dscore(ix,jx)+dx
				end if
			end do
			dscore(i,j)=rescore(i,j)
			rescore(i,j)=0.0
		end if
c
c		delete residue pair, update lenali
c
		ali1(i)=0
		ali2(j)=0
		lenali=lenali-1
		ndel=ndel+1
                ! quit bad case
                if(ndel.gt.ncand) return
	end do
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine doadditions(i,j,nnew,new,fragali1,fragali2,
     $	ali1,ali2,next1,prev1,lset,candij,ncand,d1,d2,
     $	rescore,dscore,lenali,cnt,blocksize,nacc,lsave,test,nextres,
     $	next2,nres1)
	implicit none
	include 'gagasizes.for'
	integer nnew,new(2,maxres),i,j,nres1
	integer*2 fragali1(maxres),fragali2(maxres)
        integer*2 ali1(maxres),ali2(maxres)
	integer next1(0:maxres+1),prev1(0:maxres+1),cnt(maxres)
	logical lset(maxpair)
	integer candij(2,maxpair*4),ncand,test(maxres,-maxres:maxres)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer lenali,blocksize,nacc,nextres(0:maxres+1)
        integer next2(0:maxres+1)
	logical lsave
c
	real scorefun
c
	integer inew,k,l,i1,j1,ix,jx,in,ip,iold,jold
	real s,dx
c
	if(ldebug)then
		write(*,*) 'This is doadditions',i,j,nnew,test(i,j)
		write(*,500) (new(1,k),new(2,k),k=1,nnew)
		write(*,500) (ali1(k),k=1,38)
		write(*,500) (nextres(k),k=0,39)
	end if
c
c	add gene
c
	fragali1(i)=j
	if(j.lt.0) fragali2(-j)=-i
	if(j.gt.0) fragali2(j)=i
	lset(test(i,j))=.true.
	do l=0,blocksize-1
		cnt(i+l)=cnt(i+l)+1
		if(cnt(i+l).gt.blocksize) write(*,*) 
     $			' WARNING: excess cnt',i+l,j+l
	end do
c
c	update prev1,next
c
	ip=prev1(i)
	in=next1(i)
	do ix=ip,i-1
		next1(ix)=i
	end do
	do ix=i+1,in
		prev1(ix)=i
	end do
	j1=abs(j)
	in=next2(j1)
	ix=j1
	do while((ix.gt.0).and.(next2(ix-1).eq.in))
		ix=ix-1
		next2(ix)=j1
	end do
c
c	add residues
c
	do inew=1,nnew
		i1=new(1,inew)
		j1=new(2,inew)
c
c		update nextres
c
		in=nextres(i1)
		ix=i1
c		write(*,*) 'nextres',in,i,nextres(ix-1)
		do while((ix.gt.0).and.(nextres(ix-1).eq.in))
			ix=ix-1
			nextres(ix)=i1
c			if(ldebug) write(*,*) 'set nextres',ix,' to ',i1
		end do
		if(lsave) then
c
c			update dscore(x,y), rescore(x,y)
c
			do k=1,ncand
				ix=candij(1,k)
				jx=candij(2,k)
				dx=0.0
				if((i1.ne.ix).and.(j1.ne.jx)) then
					s=scorefun(d1(ix,i1),d2(jx,j1))
					dx=s+s
				end if
				if(ali1(ix).eq.jx) then
					rescore(ix,jx)=rescore(ix,jx)+dx
				else
					dscore(ix,jx)=dscore(ix,jx)+dx
				end if
				if(ldebug) write(*,510) k,i1,j1,ix,jx,
     $					dscore(ix,jx),rescore(ix,jx),dx
			end do
			rescore(i1,j1)=dscore(i1,j1)
			dscore(i1,j1)=0.0
		end if
c
c		add residue pair, update lenali
c
		iold=ali2(j1)
		jold=ali1(i1)
		if((iold.ne.0).or.(jold.ne.0)) then
			write(*,*) ' WARNING: undeleted ',i1,j1,iold,jold
			write(*,500) (ali1(k),k=1,nres1)
			write(*,500) (fragali1(k),k=1,nres1)
		end if
		if(iold.ne.0) ali1(iold)=0
		if(jold.ne.0) ali2(jold)=0
		ali1(i1)=j1
		ali2(j1)=i1
		lenali=lenali+1
		nacc=nacc+1
	end do
c
500	format(20i4)
510	format('update',5i5,3f10.3)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine testaddition(i,j,lsave,blocksize,prev1,next1,ali1,
     $	cnt,nres1,d1,d2,rescore,dscore,dx,nrem,rem,nnew,new,ngenedel,
     $	genedel,nres2,fragali1,fragali2,nextres,next2)
	implicit none
	include 'gagasizes.for'
	integer  ngenedel,genedel(maxres)
	integer i,j,nrem,rem(2,maxres),nnew,new(2,maxres),blocksize
	integer prev1(0:maxres+1),next1(0:maxres+1),cnt(maxres)
        integer nres1,nres2
	integer*2 ali1(maxres),fragali1(maxres),fragali2(maxres)
	logical lsave
	real dx
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer nextres(0:maxres+1),next2(0:maxres+1)
c
	real getrescore
	real scorefun
c
c	in: 
c		i,j -- gene to be added
c		ltop -- TRUE means topology constraint
c	out: 
c		dx -- change in score on addition
c	    	nrem,rem -- list of residues to delete
c		nnew,new -- list of residues to add
c		ngenedel,genedel -- list of genes to delete
c
c		N.B.: nrem,nnew are initialized to zero
c
	integer k,l,inew,i1,i2
	real s
c
	if(ldebug) then
		write(*,*) 'This is testaddition',i,j,ltop,lsave
		write(*,500) (ali1(k),k=1,nres1)
	end if
c
c	initialize
c
	ngenedel=0
	nrem=0
	nnew=0
c
c	optional topology check in range <i-blocksize+1, > i+blocksize-1
c		both in i- and j-directions
c
	if(ltop)call checktopo(i,j,blocksize,prev1,next1,fragali1,
     $	ngenedel,genedel,nres1)
c
c	which genes overlap with gene (i,j) ?
c	check range i-blocksize+1 ... i+blocksize-1
c                   j-blocksize+1 ... j+blocksize-1
c
	call geneoverlap(i,j,fragali1,fragali2,next1,blocksize,
     $		nres1,nres2,ngenedel,genedel,next2)
c
c	which residues go with genedel ?
c
	call gowithgenes(i,j,ngenedel,genedel,ali1,cnt,blocksize,
     $		nrem,rem,.true.)
c
c	new residues
c
	do l=0,blocksize-1
		if(ali1(i+l).ne.abs(j+l)) then
			nnew=nnew+1
			new(1,nnew)=i+l
			new(2,nnew)=abs(j+l)
		end if
	end do
c
c	if(ldebug) then
c		write(*,*) 'genedel:',ngenedel
c		write(*,500) (genedel(k),k=1,ngenedel)
c		write(*,*) 'rem:',nrem
c		write(*,500) (rem(1,k),rem(2,k),k=1,nrem)
c		write(*,*) 'new:',nnew
c		write(*,500) (new(1,k),new(2,k),k=1,nnew)
c	end if
c
c	calculate dx
c
c	add new
c
	dx=0.0
	do inew=1,nnew
	  i1=new(1,inew)
	  i2=new(2,inew)
	  if(i1.le.maxres.and.i2.le.maxres) then
		if(lsave) then
			dx=dx+dscore(i1,i2)
                        if(ldebug) then
                          write(*,*) 'add inew',inew,i1,i2,dscore(i1,i2),dx
                        end if
		else
			dx=dx+getrescore(i1,i2,d1,d2,ali1,nextres,nres1)
		end if
c
c		add new x new
c
		do k=1,inew-1
			s=scorefun(d1(i1,new(1,k)),d2(i2,new(2,k)))
			dx=dx+s+s
		        if(ldebug) then
                         write(*,*) 'add new*new ',inew,k,i1,i2,
     $	new(1,k),new(2,k),s+s,dx,d1(i1,new(1,k)),d2(i2,new(2,k))
                        end if
		end do
c
c		subtract new x rem
c
		do k=1,nrem
			if((i1.ne.rem(1,k)).and.(i2.ne.rem(2,k))) then
				s=scorefun(d1(i1,rem(1,k)),d2(i2,rem(2,k)))
				dx=dx-s-s
                                if(ldebug) then
                           write(*,*) 'subtract new*old ',inew,k,i1,i2,
     $	new(1,k),new(2,k),-s-s,dx,d1(i1,rem(1,k)),d2(i2,rem(2,k))
                                end if
			end if
		end do
	  end if 
	end do
c
c	subtract rem x rem
c
	do k=1,nrem
		if(lsave) then
			dx=dx-rescore(rem(1,k),rem(2,k))
		else
			dx=dx-getrescore(rem(1,k),rem(2,k),d1,d2,
     $				ali1,nextres,nres1)
		end if
		do l=1,k-1
			s=scorefun(d1(rem(1,k),rem(1,l)),
     $				d2(rem(2,k),rem(2,l)))
			dx=dx+s+s
		end do
	end do
c
	if(ldebug)write(*,*) 'dx: ',dx
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine geneoverlap(i,j,fragali1,fragali2,next1,blocksize,
     $		nres1,nres2,ngenedel,genedel,next2)
	implicit none
	include 'gagasizes.for'
	integer i,j,ngenedel,genedel(maxres),blocksize,nres1,nres2,
     $		next1(0:maxres),next2(0:maxres+1)
	integer*2 fragali1(maxres),fragali2(maxres)
c
	integer k,l,j1
c
	if(ldebug)then
		write(*,*) 'This is geneoverlap',i,j
		write(*,500) (next1(k),k=1,nres1)
		write(*,500) (next2(k),k=1,nres2)
	end if
c
c	check i direction
c
	k=next1(max(0,i-blocksize))
	do while(k.le.min(nres1,i+blocksize-1))
		l=fragali1(k)
		if(k-l.ne.i-j) then
			ngenedel=ngenedel+1
			genedel(ngenedel)=k
		end if
		k=next1(k)
	end do
c
c	check j direction
c
	if(j.gt.0) then
c		parallel - parallel
		l=next2(max(0,j-blocksize))
		do while(l.le.min(nres2,j+blocksize-1))
			k=fragali2(l)
			if((k.gt.0).and.(k-l.ne.i-j).and.
     $				((k.ge.i+blocksize).or.(k.le.i-blocksize)))then
				ngenedel=ngenedel+1
				genedel(ngenedel)=k
			end if
			l=next2(l)
		end do
c		parallel - antiparallel: always clash
		l=next2(max(0,j-1))
		do while(l.le.min(nres2,j+blocksize+blocksize-1))
			k=fragali2(l)
			if((k.lt.0).and.
     $				((-k.ge.i+blocksize).or.(-k.le.i-blocksize)))then
				ngenedel=ngenedel+1
				genedel(ngenedel)=-k
			end if
			l=next2(l)
		end do
	else if(j.lt.0) then
		j1=-j
c		antiparallel - parallel: always clash
		l=next2(max(0,j1-blocksize-blocksize+1))
		do while(l.le.j1)
			k=fragali2(l)
			if((k.gt.0).and.
     $				((k.ge.i+blocksize).or.(k.le.i-blocksize)))then
				ngenedel=ngenedel+1
				genedel(ngenedel)=k
			end if
			l=next2(l)
		end do
c		antiparallel - antiparallel 
		l=next2(max(0,j1-blocksize))
		do while(l.le.min(nres2,j1+blocksize-1))
			k=fragali2(l)
			if((k.lt.0).and.(k-l.ne.i-j).and.
     $				((-k.ge.i+blocksize).or.(-k.le.i-blocksize)))then
				ngenedel=ngenedel+1
				genedel(ngenedel)=-k
			end if
			l=next2(l)
		end do
	end if
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	function getrescore(i1,i2,d1,d2,ali1,nextres,nres1)
c
c	calculates rescore(i1,i2) or dscore(i1,i2) on the fly
c
	implicit none
	include 'gagasizes.for'
	integer i1,i2,nextres(0:maxres+1),nres1
	integer*2 ali1(maxres)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real getrescore
c
	real scorefun
c
	real dx,s
	integer x,k
        integer*2 zero
c
        zero=0
	dx=scorefun(zero,zero)
	k=0
	do while(nextres(k).le.nres1)
		k=nextres(k)
c	do k=1,nres1
c	  if(ali1(k).ne.0) then
		x=ali1(k)
		if((k.ne.i1).and.(x.ne.i2)) then
			s=scorefun(d1(i1,k),d2(i2,x))
			dx=dx+s+s
		end if
c	  end if
	end do
	getrescore=dx
	if(ldebug) write(*,*) 'This was getrescore',i1,i2,dx
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine checktopo(i,j,blocksize,prev1,next1,fragali1,ngenedel,
     $		genedel,nres1)
c
c	topology check in range <i-blocksize+1, > i+blocksize-1
c				<j-blocksize+1, > j+blocksize-1
c
	implicit none
	include 'gagasizes.for'
	integer i,j,blocksize,ngenedel,genedel(maxres)
	integer prev1(0:maxres+1),next1(0:maxres+1),nres1
	integer*2 fragali1(maxres)
c
	integer k,l
c
	if(ldebug)then
		write(*,*) 'This is checktopo'
		write(*,500) ngenedel,(fragali1(k),k=1,nres1)
		write(*,500) (next1(k),k=0,nres1+1)
		write(*,500) (prev1(k),k=0,nres1+1)
	end if
	k=min(nres1+1,i+blocksize-1)
	do while(next1(k).le.nres1)
		k=next1(k)
		l=fragali1(k)
		if(l.le.j-blocksize) then
			ngenedel=ngenedel+1
			genedel(ngenedel)=k
		else
			goto 19
		end if
	end do
19	k=max(0,i-blocksize+1)
	do while(prev1(k).gt.0)
		k=prev1(k)
		l=fragali1(k)
		if(l.ge.j+blocksize) then
			ngenedel=ngenedel+1
			genedel(ngenedel)=k
		else
			goto 29
		end if
	end do
29	continue
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine testdeletion(ngenedel,genedel,
     $		blocksize,cnt,rescore,ali1,nextres,nres1,
     $		d1,d2,lsave,dx,nrem,rem)
c
c	in: 
c		ngenedel,genedel -- list of genes to be deleted
c	out: 
c		dx -- change in score on deletion
c	    	nrem,rem -- list of residues to delete
c
	implicit none
	include 'gagasizes.for'
	integer blocksize,cnt(maxres),nextres(0:maxres+1),
     $		nrem,rem(2,maxres),nres1,ngenedel,genedel(maxres)
	integer*2 ali1(maxres)
	real rescore(maxres,maxres),dx
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	logical lsave
c
	real scorefun
	real getrescore
c
	integer i,j,i1,i2,j1,j2,l,k,ix
	real s
c
c	write(*,*) 'This is testdeletion',i,j
c
c	which residues go with gene (i,j) ?
c
	do ix=1,ngenedel
	  i=genedel(ix)
	  j=ali1(i)
	  call gowithgenes(i,j,1,genedel,ali1,cnt,
     $          blocksize,nrem,rem,.false.)
c
c	  calculate dx
c
	  dx=0.0
	  do l=1,nrem
		i1=rem(1,l)
		i2=rem(2,l)
c
c		subtract rescore
c
		if(lsave) then
			dx=dx-rescore(i1,i2)
		else
			dx=dx-getrescore(i1,i2,d1,d2,ali1,nextres,nres1)
		end if
c
c		add rem x rem
c
		do k=1,l-1
			j1=rem(1,k)
			j2=rem(2,k)
			s=scorefun(d1(i1,j1),d2(i2,j2))
			dx=dx+s+s
		end do
	  end do
	end do
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine gowithgenes(i1,j1,ngenedel,genedel,ali1,cnt,blocksize,
     $		nrem,rem,laddition)
c
c	in:	ngenedel, genedel -- list of genes to delete
c		laddition -- TRUE if adding, FALSE if deleting genes
c		i1,j1 -- gene to append (don't remove identical residue pairs)
c	out:	nrem rem -- list of residues to delete (x,y always positive)
c	N.B.:	nrem is initialized to zero !
c
	implicit none
	include 'gagasizes.for'
	integer blocksize,cnt(maxres),nrem,rem(2,maxres),i1,j1
	integer ngenedel,genedel(maxres)
	integer*2 ali1(maxres)
	logical laddition
c
	integer tmp(maxres),l,k,i,m
c
c	copy cnt() to tmp()
c
	do i=1,ngenedel
		do l=0,blocksize-1
			k=abs(genedel(i)+l)
			tmp(k)=cnt(k)
		end do
	end do
c
c	subtract counts from tmp(), remove residue if empty
c
	nrem=0
	do i=1,ngenedel
		do l=0,blocksize-1
			k=abs(genedel(i)+l)
			tmp(k)=tmp(k)-1
			if(tmp(k).eq.0) then
			  if(laddition) then
c				testaddition
				m=k-i1
				if((m.ge.0).and.(m.le.blocksize-1).and.
     $					(ali1(k).eq.abs(j1+m))) goto 19
				nrem=nrem+1
				rem(1,nrem)=k
				rem(2,nrem)=ali1(k)
			  else
c				testdeletion
				nrem=nrem+1
				rem(1,nrem)=k
				rem(2,nrem)=ali1(k)
			  end if
19			end if
		end do
	end do
c
	if(ldebug) then
		write(*,*) 'This was gowithgenes',ngenedel,nrem
		write(*,500) (genedel(i),i=1,ngenedel)
		write(*,500) (rem(1,i),rem(2,i),i=1,nrem)
	end if
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine appendcand(i,j,blocksize,
     $		lset,ltest,dscore,ncand,candij,
     $		ali1,ali2,nres1,d1,d2,nextres,lsave,
     $		ntest,test,testij)
c
c	in: gene (i,j), where j may be negative
c	out: residue pair candidates (x,y), where y is always positive !
c
	implicit none
	include 'gagasizes.for'
	integer blocksize,nres1,i,j
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	logical lsave
	integer*2 ali1(maxres),ali2(maxres)
	logical lset(maxpair),ltest(maxres,maxres)
	integer ncand, candij(2,maxpair*4)
	integer ntest,testij(2,maxpair),test(maxres,-maxres:maxres),
     $		nextres(0:maxres+1)
	real dscore(maxres,maxres)
c
	real addscore
c
	integer l,j1
c
	if(ldebug) write(*,*) 'This is appendcand',i,j,ncand,lsave,ntest
c
c	add to gene-candidate list
c
	if(test(i,j).ne.0) then
		if(lverb)write(*,*) 'skip ',i,j
		return
	end if
	ntest=ntest+1
	testij(1,ntest)=i
	testij(2,ntest)=j
	test(i,j)=ntest
	lset(ntest)=.false.
	if(lsave) then
c
c	  add to residue-candidate lists
c
	  do l=0,blocksize-1
	   j1=abs(j+l)
	   if(.not.ltest(i+l,j1)) then
		ncand=ncand+1
		candij(1,ncand)=i+l
		candij(2,ncand)=j1
		ltest(i+l,j1)=.true.
c
c		calculate dscore(i+l,j1)
c
		if(ali1(i+l).ne.j1) dscore(i+l,j1)
     $			  =addscore(i+l,j1,ali1,ali2,nextres,nres1,d1,d2)
	   end if
	  end do
	end if
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	function addscore(i,j,ali1,ali2,nextres,nres1,d1,d2)
c
c	returns the marginal score of i,j relative to ali1, excluding iold,jold
c	j must be positive
c
	implicit none
	include 'gagasizes.for'
	real addscore
	integer i,j,nres1
	integer*2 ali1(maxres),ali2(maxres)
	integer nextres(0:maxres+1)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
c
	real scorefun
c
	integer*2 a,b,zero
	real s,dx
	integer k,iold,jold
c
	if(j.lt.0) write(*,*) 'This is addscore',i,j
c
        zero=0
	dx=scorefun(zero,zero)
	iold=ali2(j)
	jold=ali1(i)
	k=0
	do while(nextres(k).le.nres1)
		k=nextres(k)
		if((k.ne.iold).and.(ali1(k).ne.jold)) then
			a=d1(i,k)
			b=d2(j,ali1(k))
			s=scorefun(a,b)
			dx=dx+s+s
			if(ldebug)write(*,*) 'pair',i,j,k,ali1(k),s+s,dx
		end if
	end do
	addscore=dx

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine pullgene(gblocksize,blocksize,ngene,genepool,
     $		lset,ltest,dscore,ncand,candij,fragali1,next1,
     $		ali1,ali2,nres1,nres2,d1,d2,nextres,lsave,
     $		ntest,test,testij)
	implicit none
	include 'gagasizes.for'
	integer ngene,blocksize,gblocksize,nres1,nres2
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	logical lsave
	integer*2 genepool(5,maxpair)
	integer*2 ali1(maxres),ali2(maxres),fragali1(maxres)
	logical lset(maxpair),ltest(maxres,maxres)
	integer ncand, candij(2,maxpair*4)
	real dscore(maxres,maxres)
	integer ntest,testij(2,maxpair),test(maxres,-maxres:maxres),
     $		nextres(0:maxres+1),next1(0:maxres+1)
c
	logical lhalfin
c
	integer igene,i,j,l,p,q
c
	if(ldebug)then
		write(*,*) 'This is pullgene',ncand,ntest,nres1,nres2,ngene
		write(*,500) (ali1(i),i=1,nres1)
		write(*,500) (nextres(i),i=1,nres1)
	end if
c
	do igene=1,ngene
		do l=0,gblocksize-blocksize
			i=genepool(1,igene)+l
			j=genepool(3,igene)+l
			p=genepool(2,igene)+l
			q=genepool(4,igene)+l
c
c			gene may be antiparallel, residue pairs may not !
c
c			append one-half to cand-list if other-half is 
c				at least partly aligned in ali1
c
			if(test(p,q).eq.0) then
			  if(lhalfin(i,j,fragali1,blocksize,next1,nres1))
     $				call appendcand(p,q,blocksize,
     $				lset,ltest,dscore,ncand,candij,ali1,ali2,
     $				nres1,d1,d2,nextres,lsave,
     $				ntest,test,testij)
			end if
			if(test(i,j).eq.0) then
			  if(lhalfin(p,q,fragali1,blocksize,next1,nres1))
     $				call appendcand(i,j,blocksize,
     $				lset,ltest,dscore,ncand,candij,ali1,ali2,
     $				nres1,d1,d2,nextres,lsave,
     $				ntest,test,testij)
			end if
		end do
	end do
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	function lhalfin(i,j,fragali1,blocksize,next1,nres1)
	implicit none
	include 'gagasizes.for'
	integer i,j,blocksize,next1(0:maxres+1),nres1
	integer*2 fragali1(maxres)
	logical lhalfin
c
	integer k,l
c
	lhalfin=.false.
	k=next1(max(0,i-blocksize))
	do while((.not.lhalfin).and.(k.le.i+blocksize-1).and.(k.le.nres1))
		l=fragali1(k)
		lhalfin=(k-l.eq.i-j)
		k=next1(k)
	end do
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine printcandlist(ncand,candij,dscore,rescore)
	implicit none
	include 'gagasizes.for'
	integer ncand,candij(2,maxpair*4)
	real dscore(maxres,maxres),rescore(maxres,maxres)
c
	integer i,p,q
c
c	print candlist !
c
	write(*,*) 'This is the candidate list',ncand
	do i=1,ncand
		p=candij(1,i)
		q=candij(2,i)
		write(*,500) i,p,q,dscore(p,q),rescore(p,q)
	end do
c
500	format('candidate',3i5,2f10.3,2x,l5)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine initialize(cnt,ali1,ali2,fragali1,fragali2,
     $		prev1,next1,ltest,ncand,rescore,dscore,
     $		nres1,nres2,ntest,nextres,test,next2)
	implicit none
	include 'gagasizes.for'
	integer cnt(maxres),prev1(0:maxres+1),next1(0:maxres+1)
	integer nres1,nres2
	integer*2 ali1(maxres),ali2(maxres),fragali1(maxres)
        integer*2 fragali2(maxres)
	logical ltest(maxres,maxres)
	integer ncand, ntest,nextres(0:maxres+1),next2(0:maxres+1)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer test(maxres,-maxres:maxres)
c
	integer i,j
c
	call initali1(ali1,nres1)
	call initali1(ali2,nres2)
	call initali1(fragali1,nres1)
	call initali1(fragali2,nres2)
	do i=1,nres1
		cnt(i)=0
		do j=1,nres2
			ltest(i,j)=.false.
			rescore(i,j)=0.0
			dscore(i,j)=0.0
		end do
	end do
	do i=1,nres1
		do j=-nres2,nres2
			test(i,j)=0
		end do
	end do
	do i=0,nres1+1
		prev1(i)=0
		next1(i)=nres1+1
		nextres(i)=nres1+1
	end do
	do i=0,nres2+1
		next2(i)=nres2+1
	end do
	ncand=0
	ntest=0
c
	return
	end
c
c-------------------------------------------------------------------------------
c	
	subroutine checkcnt(cnt,fragali1,nres1,blocksize)
	implicit none
	include 'gagasizes.for'
	integer cnt(maxres),nres1,blocksize
	integer*2 fragali1(maxres)
c
	integer i,l,tmp(maxres)
	logical lp
c
	lp=.false.
	do i=1,nres1
		tmp(i)=0
	end do
	do i=1,nres1
		if(fragali1(i).ne.0) then
			do l=0,blocksize-1
				tmp(i+l)=tmp(i+l)+1
			end do
		end if
	end do
	do i=1,nres1
		if(tmp(i).ne.cnt(i)) then
			write(*,*) 'cnt error',i,cnt(i),tmp(i)
			lp=.true.
		end if
	end do
	if(lp) then
		write(*,500) (fragali1(i),i=1,nres1)
		write(*,500) (cnt(i),i=1,nres1)
	end if
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c	
