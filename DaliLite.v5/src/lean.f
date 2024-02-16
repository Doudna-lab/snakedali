c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c-------------------------------------------------------------------------------
c
	subroutine lean_mc(imax0,blocksize,nres1,nres2,d1,d2,lsave,
     $		preali1,itrimx,score,bestali1,ntetra,tetrapool)
	implicit none
	include 'gagasizes.for'
	integer imax,blocksize,nres1,nres2,itrimx
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres0)
	real totscore
	logical lsave
	integer*2 initfragali(maxres),bestali1(maxres,popsize),
     $		outfragali(maxres),tetrapool(2,maxpair)
	integer ntetra
c
c	input
c		imax0 -- max # Monte Carlo steps
c		itrimx -- trim every itrimxth step
c		blocksize -- blocksize used by alignment
c		initfragali -- set of blocksize-genes to start from
c		ltop -- TRUE means sequential constraint
c		lsave -- TRUE means use rescore,dscore; FALSE means calculate all
c	output
c		score -- scores of saved alignments
c		bestali1 -- per residue alignment
c
	integer cnt(maxres),prev1(0:maxres+1),next1(0:maxres+1)
	integer nextres(0:maxres+1),next2(0:maxres+1)
	integer*2 ali1(maxres),ali2(maxres),fragali1(maxres)
        integer*2 fragali2(maxres)
 	logical ltest(maxres,maxres)
	integer ncand, candij(2,maxpair*4)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer imax0
	real score(popsize),bestscore
c
c	local bookkeeping
c		cnt(i) -- multiplicity of genes at residue i
c		ali1(i),ali2(j) -- per residue alignment
c		fragali1(i),fragali2(j) -- per gene alignment
c		prev1(i),next1(i) -- lookup with respect to fragali1(i)
c		next2(j) -- lookup with respect to fragali2(j)
c		nextres(i) -- lookup with respect to ali1(i)
c		ltest(i,j) -- TRUE if residue pair has been appended to cand-list
c		ncand -- length of candidate (residue pairs i,j) list
c		candij(2,ix) -- pointer from candidate to residue pair
c		rescore(i,j) -- worth of residue pair in current alignment
c		dscore(i,j) -- how much residue pair would add to current ali
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
	real alfa,delta
	real dx,p
	integer i,j,ix
	integer nrem,rem(2,maxres),nnew,new(2,maxres)
	integer di(maxpair),ve(maxpair),q,imp,ngenedel,genedel(maxres)
	real impscore
c
	integer jnk1(dim1)
	common /junk1/jnk1
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23(1),rescore)
	equivalence (jnk23(maxres*maxres+1),dscore)
	integer jnk5(dim5)
	common /junk5/jnk5
	equivalence (jnk5(1),candij)
	equivalence (jnk5(maxpair*8+1),di)
	equivalence (jnk5(maxpair*9+1),ve)
	integer jnk6(dim6)
	common /junk6/jnk6
	equivalence (jnk6,ltest)
c
	imax=imax0
	if(lverb)write(*,*) 'This is lean Monte Carlo',imax,ntetra,
     $		ltop,lsave,nres1,nres2
c
c	initialize
c
!	totscore=bestscore
!	call gettotscore(preali1,d1,d2,nres1,totscore)
!	if(totscore.gt.bestscore) then
!		do i=1,nres1
!			bestali1(i,1)=preali1(i)
!		end do
!		bestscore=totscore
!	end if

	bestscore=0.0
	call initialize_lean(cnt,ali1,ali2,fragali1,fragali2,prev1,
     $		next1,ltest,ncand,rescore,dscore,nres1,nres2,nextres,
     $		next2)
c
c	load candidates
c
	do ix=1,ntetra
		i=tetrapool(1,ix)
		j=tetrapool(2,ix)
		call appendcand_lean(i,j,blocksize,ltest,dscore,
     $	ncand,candij,ali1,ali2,nres1,d1,d2,nextres,lsave)
	end do
	write(*,*) 'ncand:',ncand
	write(*,500) (preali1(i),i=1,min(maxres,nres1))
c
c	convert preali1 to initfragali -- !!! assumes ltop=.TRUE. !!!
c
	do i=1,min(maxres,nres1)
		initfragali(i)=0
	end do
	do i=1,min(maxres,nres1)-3
		j=preali1(i)
		if(j.gt.0.and.j.le.min(maxres,nres2)-3) then
c
c			at least one candidate per segment, more on long ones
c
			if(i.gt.1) then
				if(preali1(i-1).eq.0) initfragali(i)=j
			end if
			if(preali1(i+3).eq.j+3) initfragali(i)=j
		end if
	end do
c
c	load initali
c
	lenali=0
	totscore=0.0
	do i=1,min(maxres,nres1)
		j=initfragali(i)
		if(j.ne.0) then
		  call testaddition(i,j,lsave,blocksize,prev1,next1,
     $			ali1,cnt,nres1,d1,d2,rescore,dscore,
     $			dx,nrem,rem,nnew,new,ngenedel,genedel,
     $			nres2,fragali1,fragali2,nextres,next2)
		  ! aconitase fix 20-Sep-95				    !
		  ! force accept even if dx=0 from dscore < short cand-list !
		  bestscore=0.0
		  call changeconfig_lean(ngenedel,genedel,nrem,rem,i,j,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)
		  !write(*,*),i,j,nrem,nnew,bestscore,totscore
		end if
	end do
	if(lverb)write(*,*) 'initial totscore',totscore,lenali
cwrite(*,500) (bestali1(i,1),i=1,nres1)
!	call gettotscore(ali1,d1,d2,nres1,totscore)
!        write(*,*) 'gettotscore',totscore
c	write(*,500) (initfragali(i),i=1,nres1)
!	write(*,500) (ali1(i),i=1,nres1)
c
c	Monte Carlo step
c
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
c		adding phase  WARM
c
c
c		test in randomized order
c
		do ix=1,ntetra
			di(ix)=ix
			ve(ix)=nint(100.0*ran(seed))
		end do
		call j_index_Qsort(1,ntetra,maxpair,di,ve)
		do q=1,ntetra
			ix=di(q)
			i=tetrapool(1,ix)
			j=tetrapool(2,ix)
			if(fragali1(i).ne.j) then
				call testaddition(i,j,lsave,blocksize,
     $	prev1,next1,ali1,cnt,nres1,d1,d2,
     $	rescore,dscore,dx,nrem,rem,nnew,new,
     $	ngenedel,genedel,nres2,
     $	fragali1,fragali2,nextres,next2)
				p=getp(alfa*dx)
				if(ran(seed).lt.p) then
c				  write(*,*)'accepted',i,j,dx

				  call changeconfig_lean(ngenedel,
     $	genedel,nrem,rem,i,j,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)
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
!		if(mod(istep,itrimx).eq.0) then
                if(istep.gt.0) then
		 i=0
		 do while(next1(i).le.nres1)
		  i=next1(i)
		  j=fragali1(i)
c		  Trim from ends only
		  if((cnt(i).eq.1).or.(cnt(i+blocksize-1).eq.1))then
			ngenedel=1
			genedel(1)=i
			call testdeletion(ngenedel,genedel,blocksize,
     $	cnt,rescore,ali1,nextres,nres1,d1,d2,lsave,dx,
     $	nrem,rem)
c			p=getp(alfa*dx)
c			if(ran(seed).lt.p) then
c			only accept improvements !
			if(dx.gt.0.0) then
			  call changeconfig_lean(ngenedel,
     $	genedel,nrem,rem,0,0,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,blocksize,lsave,totscore,bestscore,
     $	nres1,outfragali,bestali1,dx)
				  if(totscore.gt.impscore+1.0) then
					imp=istep
					impscore=totscore
				  end if

			end if
		  end if
		 end do
		end if
		if(lverb)write(*,510) istep,imp,nacc,ndel,lenali,
     $			ncand,totscore,bestscore,alfa
                score(1)=bestscore/100 ! case return eearly
                ! quit hopelss cases
                if(nacc.gt.ncand) return
                if(ndel.gt.ncand) return
c!!!		call gettotscore(ali1,d1,d2,nres1,totscore)
c		if(lsave)call print(ncand,candij,dscore,rescore,ali1,
c     $			next1,d1,d2,nres1)
		if(istep-imp.ge.20) istep=imax0+1
		if(nacc/float(ntetra).gt.0.1) alfa=10*alfa
	end do
	score(1)=bestscore/100
	write(*,*) 'lean bestali1 ',bestscore
	write(*,500) (bestali1(i,1),i=1,nres1)
c
500	format(20i4)
510	format('step',6i6,2f10.2,f10.5)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine changeconfig_lean(ngenedel,genedel,nrem,rem,i,j,
     $	nnew,new,fragali1,fragali2,ali1,ali2,next1,prev1,next2,
     $	nextres,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,ndel,nacc,blocksize,lsave,totscore,bestscore,
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
	integer cnt(maxres),ndel,nacc,blocksize
	logical lsave
	integer nres1
	integer*2 outfragali(maxres),bestali1(maxres)
c
	integer k
c
	call dodeletions_lean(ngenedel,genedel,nrem,rem,fragali1,fragali2,
     $	ali1,ali2,next1,prev1,candij,ncand,d1,d2,rescore,dscore,
     $	lenali,cnt,ndel,lsave, nextres,blocksize,next2)
	if(i.ne.0) call doadditions_lean(i,j,nnew,new,fragali1,fragali2,
     $	ali1,ali2,next1,prev1,candij,ncand,d1,d2,rescore,dscore,lenali,
     $	cnt,blocksize,nacc,lsave,nextres,next2,nres1)
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
		call saveali(nres1,fragali1,outfragali,totscore,
     $                  bestscore)
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
	subroutine dodeletions_lean(ngenedel,genedel,nrem,rem,fragali1,
     $		fragali2,ali1,ali2,next1,prev1,candij,ncand,d1,d2,
     $		rescore,dscore,lenali,cnt,ndel,lsave,nextres,blocksize,
     $		next2)
	implicit none
	include 'gagasizes.for'
	integer ngenedel,genedel(maxres),nrem,rem(2,maxres)
	integer next2(0:maxres+1)
	integer*2 fragali1(maxres),fragali2(maxres),ali1(maxres)
        integer*2 ali2(maxres)
	integer next1(0:maxres+1),prev1(0:maxres+1),cnt(maxres)
	integer candij(2,maxpair*4),ncand,blocksize
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real rescore(maxres,maxres),dscore(maxres,maxres)
	integer lenali,ndel,nextres(0:maxres+1)
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
                        if(ix.le.0) exit
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
                        if(ix.le.0) exit
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
	end do
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine doadditions_lean(i,j,nnew,new,fragali1,fragali2,
     $		ali1,ali2,next1,prev1,candij,ncand,d1,d2,
     $		rescore,dscore,lenali,cnt,blocksize,nacc,lsave,nextres,
     $		next2,nres1)
	implicit none
	include 'gagasizes.for'
	integer nnew,new(2,maxres),i,j,nres1
	integer*2 fragali1(maxres),fragali2(maxres),ali1(maxres)
        integer*2 ali2(maxres)
	integer next1(0:maxres+1),prev1(0:maxres+1),cnt(maxres)
	integer candij(2,maxpair*4),ncand
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
		write(*,*) 'This is doadditions',i,j,nnew
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
        if(ix.gt.0) then
	  do while(next2(ix-1).eq.in)
		ix=ix-1
		next2(ix)=j1
                if(ix.le.0) exit
	  end do
        end if
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
                        if(ix.le.0) exit
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
	subroutine appendcand_lean(i,j,blocksize,
     $		ltest,dscore,ncand,candij,
     $		ali1,ali2,nres1,d1,d2,nextres,lsave)
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
	logical ltest(maxres,maxres)
	integer ncand, candij(2,maxpair*4)
	integer nextres(0:maxres+1)
	real dscore(maxres,maxres)
c
	real addscore
c
	integer l,j1
c
	if(ldebug) write(*,*) 'This is appendcand',i,j,ncand,lsave
	if(lsave) then
	  do l=0,blocksize-1
	   j1=abs(j+l)
	   if(.not.ltest(i+l,j1)) then
c
c		add to residue-candidate list
c			used to update rescore(x,y),dscore(x,y)
c
		ncand=ncand+1
		candij(1,ncand)=i+l
		candij(2,ncand)=j1
		ltest(i+l,j1)=.true.
c
c		calculate dscore(i+l,j1)
c
		if(ali1(i+l).ne.j1) dscore(i+l,j1)
     $		  =addscore(i+l,j1,ali1,ali2,nextres,nres1,d1,d2)
	   end if
	  end do
	end if
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine initialize_lean(cnt,ali1,ali2,fragali1,fragali2,
     $		prev1,next1,ltest,ncand,rescore,dscore,
     $		nres1,nres2,nextres,next2)
	implicit none
	include 'gagasizes.for'
	integer cnt(maxres),prev1(0:maxres+1),next1(0:maxres+1)
	integer nres1,nres2
	integer*2 ali1(maxres),ali2(maxres),fragali1(maxres)
        integer*2 fragali2(maxres)
	logical ltest(maxres,maxres)
	integer ncand, nextres(0:maxres+1),next2(0:maxres+1)
	real rescore(maxres,maxres),dscore(maxres,maxres)
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
	do i=0,nres1+1
		prev1(i)=0
		next1(i)=nres1+1
		nextres(i)=nres1+1
	end do
	do i=0,nres2+1
		next2(i)=nres2+1
	end do
	ncand=0
c
	return
	end
c
c-------------------------------------------------------------------------------
c
