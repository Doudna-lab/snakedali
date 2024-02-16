c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
	subroutine getalignments(code1,chainid1,code2,chainid2,d1,d2,
     $		nres1,nres2,gblocksize,blocksize,ngene,genepool,
     $		ali1,score,lpreali,preali1)
        implicit none
	include 'gagasizes.for'
	character chainid1,chainid2
	character*4 code1,code2
	integer nres1,nres2,blocksize,gblocksize
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
	integer*2 ali1(maxres,popsize)
	real score(popsize)
	logical lpreali
c
	character*80 pairfilnam
c
	character*80 filnam
c
c	genepool: 1-4: i1 i2 j1 j2 5: score 
c
	integer*2 genepool(5,maxpair)
	integer ngene,i
c
	filnam=pairfilnam(code1,chainid1,code2,chainid2,'ali ')
	if(lverb)write(*,*) filnam,nres1,nres2
	call simple(ngene,genepool,nres1,nres2,gblocksize,blocksize,
     $		score,ali1,d1,d2,lpreali,preali1)
	if(lverb)write(*,*) score(1)
	if(lverb)write(*,500) (ali1(i,1),i=1,nres1)
c
500	format(20i4)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	function hack(n)
	implicit none
	include 'gagasizes.for'
	real hack
	integer n
c
c	define minimal score required to do refinement
c	e.g. min(nres,100) if d0=0.20
c
	hack=min(hack0,hack1*n)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine simple(ngene,genepool,nres1,nres2,gblocksize,
     $		blocksize,score,bestfragali1,d1,d2,lpreali,preali1)
	implicit none
	include 'gagasizes.for'
	integer ngene,nres1,nres2,gblocksize,blocksize
	integer*2 genepool(5,maxpair)
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
	real score(popsize)
	integer*2 bestfragali1(maxres,popsize),bestali1(maxres,popsize)
	logical lpreali
c
	real ran, hack
c
	integer ifr
	real bestscore0(maxpair),totscore0
	integer*2 initali(maxres),initfragali(maxres,popsize),
     $		outfragali(maxres)
c
	integer itrip
	logical lfrozen(maxpair),lf,lkill
c
	integer i,ix,j,k,l
	real topscore,scorecutoff,round,expt,tmpscore(popsize)
	integer nblock,block(maxres)
	logical lkeep(maxres),lb,lempty
	integer di(popsize)
c
c	ldebug=.false.
c
c	build-up all seeds
c
	call buildup(ngene,genepool,nres1,nres2,blocksize,score,lfrozen,
     $		d1,d2,bestfragali1,lempty,lpreali,preali1)
	if(lempty) then
		if(lverb)write(*,*) 'no seeds - skip alignment'
		return
	end if
	topscore=0.0
	do i=1,popsize
		bestscore0(i)=score(i)*100
		do ix=1,nres1
			bestali1(ix,i)=0
		end do
	end do
	call sortscores(di,score)
	write(*,*) ' sorted scores'
	do i=1,popsize
		write(*,*) i,di(i),score(di(i))
	end do
	topscore=score(di(nbest))
c
c	now update bestfragali1 ! to account for change geneblocksize->blocksize
c
	if(lverb)write(*,*) 'geneblocksize:',gblocksize,
     $		' blocksize:',blocksize
	do itrip=1,popsize
		call initali1(initfragali,nres1)
		do i=1,nres1
			j=bestfragali1(i,itrip)
			if(j.ne.0) then
				do l=0,gblocksize-blocksize
					initfragali(i+l,1)=j+l
				end do
			end if
		end do
		do i=1,nres1
			bestfragali1(i,itrip)=initfragali(i,1)
		end do
	end do
c
c     optimize popsize best seeds
c
      scorecutoff=0.0
      round=1.0
      ifr=1
      lf=.false.
      do while(.not.lf)
	lf=.true.
c	freeze alignments with scores lower than scorecutoff
	ifr=ifr+20
	round=round+1.0
	expt=100.0*nres1*(1.0-1.0/round)
	scorecutoff=(1.0-1.0/round)*topscore
	do itrip=1,popsize
	  if((score(itrip).gt.0.0).and.(score(itrip).lt.scorecutoff)) then
	    if(lverb)write(*,*)' too low score: killing ',
     $		itrip,score(itrip),scorecutoff
	    lfrozen(itrip)=.true.
	    score(itrip)=0.0
	  end if
	end do
	do itrip=1,popsize
	 if(.not.lfrozen(itrip)) then
	  if(lverb)write(*,*) ' seed',itrip,lfrozen(itrip),
     $          score(itrip),ifr
	  do i=1,nres1
		initali(i)=0
		initfragali(i,1)=bestfragali1(i,itrip)
	  end do
	  i=20
	  call mclean(i,ngene,genepool,gblocksize,blocksize,nres1,nres2,
     $		d1,d2,totscore0,.true.,initfragali,initali,itrim(2),ikill(2),
     $		bestali1,score,killcutoff(2),itrip,lkill,lfrozen,
     $		outfragali,0,nclone(2),expt)
	  lf=.false.
	  if(ldebug)write(*,500) (initali(i),i=1,nres1)
c!!!	  call gettotscore(initali,d0,d1,d2,nres1,ltop,totscore0)
c
c	  kill lower-score converged alignments
c	  save alignment in bestali1,bestfragali1
c
	  if(totscore0.gt.bestscore0(itrip)+1.0) then
		  if(.not.lfrozen(itrip)) then
			do i=1,nres1
				bestfragali1(i,itrip)=outfragali(i)
				bestali1(i,itrip)=initali(i)
			end do
			score(itrip)=totscore0/100
			bestscore0(itrip)=totscore0
		  end if
	  else
		lfrozen(itrip)=.true.
	  end if
	  if(lkill) then
		score(itrip)=0.0
		lfrozen(itrip)=.true.
		if(lverb)write(*,*) ' killing converged', itrip
 	  end if
c
c!	  require at least expected score
c
	  if(totscore0.lt.expt) then
		lfrozen(itrip)=.true.
		if(lverb)write(*,*) 'too low score',totscore0,' killing ',itrip
		totscore0=0.0
	  end if
	  write(*,*) ' cycle ',itrip,score(itrip),lfrozen(itrip)
	 end if
	 call sortscores(di,score)
	 topscore=score(di(nbest))
	 lf=(lf.and.lfrozen(itrip))
	end do
c	end of seed
	write(*,*) ' end of cycle ',lf
c
      end do
c     until all seeds frozen
	call sortscores(di,score)
	write(*,*) ' sorted scores'
	do i=1,popsize
		write(*,*) i,di(i),score(di(i))
	end do
	write(*,*) ' best cycle was ',di(1),score(di(1))
c
c	refine Nbest alignments
c
	do ix=1,nbest
	  write(*,*) ' refine ',ix,di(ix),score(di(ix))
c
c	  a bit of a hack
c
	  if(score(di(ix)).lt.hack(nres1)) then
	    if(lverb)write(*,*) ' too low score -- skip refinement'
	    goto 199
	  end if
c
c	  while score improves
c
	  lf=.false.
	  topscore=score(di(ix))*100
c
c	  initialize
c
	  do i=1,popsize
		bestscore0(i)=0.0
		do k=1,nres1
			bestali1(k,i)=0
		end do
	  end do
	do while(.not.lf)
c
c	 generate popsize initial alignments that keep 70 % of topfragali blocks
c	 initialize bestali1 for killer
c
	 call getblocks(bestfragali1,di(ix),blocksize,block,nblock,nres1)
	 if(nblock.le.1) then
		write(*,*) nblock,' is too few blocks for refinement '
		goto 199
	 end if
	 l=0
	 do i=1,nres1
		if(bestali1(i,1).ne.0) l=l+1
		initfragali(i,1)=bestfragali1(i,di(ix))
	 end do
         do itrip=2,popsize
		lfrozen(itrip)=.false.
c		must select at least one block
		lb=.false.
		do while(.not.lb)
		  do i=1,nblock
			lkeep(i)=(ran(seed).lt.pblock)
			lb=(lb.or.lkeep(i))
		  end do
		end do
		do i=1,nres1
			initfragali(i,itrip)=0
			k=block(i)
			if(k.gt.0) then
			 if(lkeep(k)) initfragali(i,itrip)=bestfragali1(i,di(ix))
			end if
		end do
	 end do
	 do itrip=1,popsize
		tmpscore(itrip)=0.0
	 end do
	 lf=.true.
	 i=20
	 call mclean(i,ngene,genepool,gblocksize,blocksize,nres1,nres2,
     $		d1,d2,totscore0,.true.,initfragali,initali,itrim(3),ikill(3),
     $		bestali1,tmpscore,killcutoff(3),1,lkill,lfrozen,
     $		outfragali,l,nclone(3),0.0)
	 if(lverb)write(*,500) (initali(i),i=1,nres1)
	 if(totscore0.ge.topscore+refitol) then
		  if(lverb)write(*,*) 'save alignment',ix,di(ix)
		  do i=1,nres1
			bestfragali1(i,di(ix))=outfragali(i)
			bestali1(i,1)=initali(i)
		  end do
		  score(di(ix))=totscore0/100
		  topscore=totscore0
		  lf=.false.
	  end if
	  write(*,*) ' score after melt',totscore0,topscore,lf,ix,
     $          score(di(ix))
c!!!	  call gettotscore(initali,d1,d2,nres1,totscore0)
	end do
c	end of Nbest refinements
199	end do
c
	call sortscores(di,score)
	write(*,*) 'best alignment',di(1),score(di(1))
	if(lverb)write(*,500) (bestali1(i,di(1)),i=1,nres1)
	if(lverb)write(*,500) (bestfragali1(i,di(1)),i=1,nres1)
c
500	format(20i4)
c
	return
	end
c
c       ---------------------------------------------------------------------
c
	subroutine buildup(ngene,genepool,nres1,nres2,blocksize,score,
     $		lfrozen,d1,d2,bestfragali1,lempty,lpreali,preali1)
c
c	bestfragali1 returns per gene alignment
c
        implicit none
	include 'gagasizes.for'
	integer ngene,nres1,nres2,blocksize
	integer*2 genepool(5,maxpair)
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
	real score(popsize)
	integer*2 bestali1(maxres,popsize),bestfragali1(maxres,popsize)
	logical lfrozen(maxpair),lempty,lpreali
c
	integer findplace
c
	integer ntrip,trip(3,10000),singletseed(3,maxpair),ns
	real totscore
	integer*2 ali1(maxres)
	integer i,j,l
	real totscore0
	integer*2 initali(maxres),outfragali(maxres)
	integer itrip
	logical lkill
	integer ix,ifr
c
	logical ltest
c
	integer jnk4(dim4)
	common /junk4/jnk4
	equivalence (jnk4(1),trip)
	equivalence (jnk4(30001),singletseed)
c
	if(lpreali) then
	  ntrip=1
	  call initali1(ali1,nres1)
	  do i=1,nres1
		initali(i)=preali1(i)
	  end do
	  do i=1,nres1
	    if((initali(i).gt.0).and.(initali(i+3).eq.initali(i)+3)) then
		lempty=.false.
		ali1(i)=initali(i)
	    end if
	  end do
	else
	  call triplet(ngene,genepool,trip,ns)
	  if(lverb)write(*,*) ns,' triplets'
	  call tripletcontigs(nres1,nres2,trip,ns,genepool,singletseed,
     $		ntrip,blocksize)
	  if(lverb)write(*,*) ntrip, ' seeds'
	  do i=1,popsize
		score(i)=0.0
		do j=1,nres1
			bestfragali1(j,i)=0
		end do
	  end do
	end if
	lempty=(ntrip.eq.0)
	if(lempty) return
c
c!!!testing
c
	ltest=.false.
	if(ltest) then
		ntrip=1
		write(*,*) 'enter seed: i,j,len'
		read(*,*) i,j,l
		call initali1(initali,nres1)
		call initali1(ali1,nres1)
		do ix=0,l-blocksize
			ali1(i+l)=j+l
		end do
	end if
	do itrip=1,ntrip
		lfrozen(itrip)=.false.
		if(lpreali) then
			write(*,*) 'seed from prealigment'
		else if(.not.ltest) then
		    	call init2(itrip,singletseed,ali1,nres1,blocksize)
			call initali1(initali,nres1)
			do i=1,nres1
				if(ali1(i).ne.0) then
					do ix=0,blocksize-1
						initali(i+ix)=ali1(i)+ix
					end do
				end if
			end do
		end if
		totscore=0.0
		if(ldebug)then
			write(*,500) (initali(ix),ix=1,nres1)
	  		call gettotscore(initali,d1,d2,nres1,totscore)
		end if
		lkill=.false.
		if(itrip.gt.popsize) call 
     $			killer(bestali1,initali,totscore/100,nres1,blocksize,
     $				score,killcutoff(1),lfrozen,lkill,itrip,0)
		if(lkill) then
			lfrozen(itrip)=.true.
			score(itrip)=0.0
			if(lverb)write(*,*) ' weeding converged', itrip
		else
		  ifr=1
		  call mclean(ifr,ngene,genepool,blocksize,blocksize,nres1,
     $			nres2,d1,d2,totscore,.false.,ali1,initali,itrim(1),
     $			ikill(1),bestali1,score,killcutoff(1),itrip,
     $			lkill,lfrozen,outfragali,0,nclone(1),0.0)
		  if(lverb)write(*,*) ' startup:',itrip,totscore
		end if
c	  	call gettotscore(initali,d1,d2,nres1,totscore)
		totscore0=totscore
		if(lkill) then
			lfrozen(itrip)=.true.
			score(itrip)=0.0
			if(lverb)write(*,*) ' killing converged', itrip
		end if
		ix=itrip
		if(.not.lfrozen(itrip)) then
		  if(itrip.gt.popsize) then
			ix=findplace(totscore0/100,score)
			if(ix.eq.0) then
			  if(lverb)write(*,*) 'too low score, killing',itrip
			  lfrozen(itrip)=.true.
			else
			  if(lverb)write(*,*) ' overwriting ',ix,' by ',itrip
			  lfrozen(ix)=lfrozen(itrip)
			end if
		  end if
		end if
c
c		save per-residue alignment in bestali1 - for killer
c		save per-hexa alignment in bestfragali1 
c
		if(.not.lfrozen(itrip)) then
		  do i=1,nres1
			bestali1(i,ix)=initali(i)
			bestfragali1(i,ix)=outfragali(i)
		  end do
		  score(ix)=totscore0/100
		end if
	end do
c
	do itrip=1,popsize
	  if(lverb)write(*,*) 'buildup:',itrip,score(itrip),lfrozen(itrip)
	end do
c
500	format(20i4)
c
	return
	end
c
c       ---------------------------------------------------------------------
c
	function findplace(testscore,score)
c
c	returns 0 if testscore<all score(i) or ix if testscore>score(ix)
c
	implicit none
	include 'gagasizes.for'
	integer findplace
	real testscore,score(popsize)
c
	integer ix,i
	real smin
c
	ix=1
	smin=score(ix)
	do i=2,popsize
		if(score(i).lt.smin) then
			smin=score(i)
			ix=i
		end if
	end do
	if(testscore.lt.smin) then
		findplace=0
	else
		findplace=ix
	end if
c
	return
	end
c
c       ---------------------------------------------------------------------
	subroutine saveali(nres1,ali1,initali,totscore,totscore0)
	implicit none
	include 'gagasizes.for'
	integer nres1
	real totscore,totscore0
	integer*2 ali1(maxres),initali(maxres)
	integer i
c
	do i=1,nres1
		initali(i)=ali1(i)
	end do
	totscore0=totscore	
c
	return
	end
c
c-------------------------------------------------------------------------------
c
        subroutine killer(ali1,initali,testscore,
     $		nres1,blocksize,score,simcut,lfrozen,lkill,ind,reflen)
        implicit none
        include 'gagasizes.for'
        integer nres1,blocksize,kill(popsize),nkill,ind,reflen
        real score(popsize),simcut,testscore
        integer*2 ali1(maxres,popsize),initali(maxres)
	logical lfrozen(maxpair),lkill
c
        integer n,i,j,i1
        real sim(popsize)
c
c       kill lower-score alis which are more than simcut similar
c       similarity === shared residues / length of ind-ali
c	lkill=.true. means kill newcomer
c
	lkill=.false.
        nkill=0
        do i=1,popsize
                sim(i)=0.0
        end do
        n=0
        do j=1,nres1
          i1=initali(j)
          if(i1.ne.0) then
                n=n+1
                do i=1,popsize
		  if((.not.lfrozen(i)).and.(i.ne.ind)) then
                        if(ali1(j,i).eq.i1) sim(i)=sim(i)+1.0
		  end if
                end do
          end if
        end do
c
c	set n to longer of ind-lenali or reference lenali (reflen)
c
 	n=max(n,reflen)
        if(n.eq.0) return
        do i=1,popsize
c	if(.not.lfrozen(i))write(*,*) 'killer:',i,sim(i),n,score(i),testscore,
c     $		lkill,ind,sim(i)/n,simcut,reflen
          if((sim(i)/n.ge.simcut).and.(score(i).gt.0.0).and.
     $		(i.ne.ind))then
                if(score(i).lt.testscore) then
			j=i
                	nkill=nkill+1
                	kill(nkill)=j
		else
			lkill=.true.
			return
		end if
          end if
        end do
        do j=1,nkill
           if(lverb)write(*,*) ' killing ',kill(j),' by ',ind,
     $		score(kill(j)),testscore
           score(kill(j))=0.0
	   lfrozen(kill(j))=.true.
        end do
 
        return
        end
c
c-------------------------------------------------------------------------------
	subroutine init2(itrip,singletseed,initali,nres1,blocksize)
	implicit none
	include 'gagasizes.for'
	integer singletseed(3,maxpair),nres1,itrip,blocksize
	integer*2 initali(maxres)
c
	integer ix,i
c
	    ix=itrip
	    do i=1,nres1
		initali(i)=0
	    end do
	    do i=0,singletseed(3,ix)-blocksize
		initali(singletseed(1,ix)+i)=singletseed(2,ix)+i
	    end do
	    if(lverb)write(*,*) ' seed: ',ix,(singletseed(i,ix),i=1,3)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine initali1(ali1,nres)
	integer nres
	integer*2 ali1(nres)
c
	integer i
c
	do i=1,nres
		ali1(i)=0
	end do
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine tripletcontigs(nres1,nres2,trip,ntrip,genepool,
     $		singletseed,n,blocksize)
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2,ntrip,trip(3,10000),blocksize
	integer*2 genepool(5,maxpair)
c
	integer i,j,ix,k,l,n,l1,l2
	real map(maxres,-maxres:maxres),smax
c
	integer singletseed(3,maxpair)
c
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23,map)
c
	if(lverb)write(*,*) 'This is tripletcontigs',
     $		nres1,nres2,ntrip,n,blocksize
c
c	postplot cumulative triplet scores
c
	do i=1,nres1
		do j=-nres2,nres2
			map(i,j)=0.0
		end do
	end do
	do i=1,ntrip
	  do ix=1,3
		j=trip(ix,i)
		k=genepool(1,j)
		l=genepool(3,j)
		map(k,l)=map(k,l)+genepool(5,j)
		if(map(k,l).gt.smax) smax=map(k,l)
		j=trip(ix,i)
		k=genepool(2,j)
		l=genepool(4,j)
		map(k,l)=map(k,l)+genepool(5,j)
		if(map(k,l).gt.smax) smax=map(k,l)
	  end do
	end do
	if(lverb)write(*,*) ' max tripletscore ',smax
	smax=smax*seedcutoff
	n=0
	do i=1,nres1
		do j=-nres2,nres2
			if(map(i,j).le.smax) then
				map(i,j)=0.0
				n=n+1
			end if
		end do
	end do
	if(lverb)write(*,*) n,' low-score pairs blanked out'
c	if(lplottriplets) call plottriplets(nres1,nres2)
c
c	make list of contigs along each diagonal
c
	n=0
c	positive half
	do i=2,nres1
		k=i-1
		l1=0
		l2=0
		do j=1,nres2
			if(j+k.le.nres1) then
				if(map(j+k,j).gt.0.0) then
					if(l1.eq.0) l1=j+k
					l2=j+k+blocksize-1
					if(l2.eq.nres1) then
						n=n+1
	if(lverb)write(*,500) n,l1,l2,l1-k,l2-k
	singletseed(1,n)=l1
	singletseed(2,n)=l1-k
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				else
					if((j+k.gt.l2).and.(l1.gt.0)) then
						n=n+1
	if(lverb)write(*,500) n,l1,l2,l1-k,l2-k
	singletseed(1,n)=l1
	singletseed(2,n)=l1-k
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				end if
			end if
		end do
	end do
	do j=1,nres2
		k=1-j
		l1=0
		l2=0
		do i=1,nres1
			if(i-k.le.nres2) then
				if(map(i,i-k).gt.0.0) then
					if(l1.eq.0) l1=i-k
					l2=i-k+blocksize-1
					if(l2.eq.nres2) then
						n=n+1
	if(lverb)write(*,500) n,l1+k,l2+k,l1,l2
	singletseed(1,n)=l1+k
	singletseed(2,n)=l1
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				else
					if((i-k.gt.l2).and.(l1.gt.0)) then
						n=n+1
	if(lverb)write(*,500) n,l1+k,l2+k,l1,l2
	singletseed(1,n)=l1+k
	singletseed(2,n)=l1
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				end if
			end if
		end do
	end do
c	negative half
	do i=2,nres1
		k=i+nres2
		l1=0
		l2=0
		do j=-nres2,-blocksize
			if(j+k.le.nres1) then
				if(map(j+k,j).gt.0.0) then
					if(l1.eq.0) l1=j+k
					l2=j+k+blocksize-1
					if(l2.eq.nres1) then
						n=n+1
	if(lverb)write(*,500) n,l1,l2,l1-k,l2-k
	singletseed(1,n)=l1
	singletseed(2,n)=l1-k
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				else
					if((j+k.gt.l2).and.(l1.gt.0)) then
						n=n+1
	if(lverb)write(*,500) n,l1,l2,l1-k,l2-k
	singletseed(1,n)=l1
	singletseed(2,n)=l1-k
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				end if
			end if
		end do
	end do
	do j=-nres2,-blocksize
		k=1-j
		l1=0
		l2=0
		do i=1,nres1
			if(i-k.le.-1) then
				if(map(i,i-k).gt.0.0) then
					if(l1.eq.0) l1=i-k
					l2=i-k+blocksize-1
					if(l2.eq.nres2) then
						n=n+1
	if(lverb)write(*,500) n,l1+k,l2+k,l1,l2
	singletseed(1,n)=l1+k
	singletseed(2,n)=l1
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				else
					if((i-k.gt.l2).and.(l1.gt.0)) then
						n=n+1
	if(lverb)write(*,500) n,l1+k,l2+k,l1,l2
	singletseed(1,n)=l1+k
	singletseed(2,n)=l1
	singletseed(3,n)=l2-l1+1
						l1=0
						l2=0
					end if
				end if
			end if
		end do
	end do
c
500	format(' contig ',i5,2(i5,' -',i5))	
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine plottriplets(nres1,nres2)
	implicit none
	include 'gagasizes.for'
	character dummy(maxres)
	real map(maxres,-maxres:maxres),psmap(maxres,maxres)
	integer i,j,k,nres1,nres2
c
	integer jnk23(dim23)
	common /junk23/jnk23
	equivalence (jnk23,map)
c
	do i=1,nres1
		dummy(i)=' '
		do j=1,nres2
			if(map(i,j).gt.0) then
				do k=0,geneblocksize-1
					psmap(i+k,j+k)=1
				end do
			end if
		end do
	end do
c	call postplot(psmap,nres1,nres2,.true.,.false.,dummy,maxres,95,96)
c
	return
	end
c
c-------------------------------------------------------------------------------
c
        subroutine gettotscore(ali1,d1,d2,nres1,xtotscore)
        implicit none
        include 'gagasizes.for'
	real totscore,xtotscore
        integer nres1
        integer*2 d1(maxres,maxres),d2(maxres,maxres)
        integer*2 ali1(maxres)
c
        real scorefun
c
        integer i,j,k,l,q,r,a(maxres),n
        real s(maxres),x
c
        totscore=0.0
        do i=1,nres1
                s(i)=0.0
        end do
	n=0
	do i=1,nres1
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
			x=scorefun(d1(k,l),d2(q,r))
			s(i)=s(i)+x
			if(i.ne.j) s(j)=s(j)+x
			totscore=totscore+x
		end do
	end do
        if(ldebug) then
                do i=1,nres1
                  if(ali1(i).ne.0) write(*,*) 'score ',i,s(i)
                end do
        end if
        if(abs(totscore-xtotscore).gt.1.0) then
              if(xtotscore.ne.0.0) then
		write(*,*) ' !!! BUG !!! resetting totscore ',xtotscore,
     $			totscore,xtotscore-totscore,n
		write(*,500) (ali1(i),i=1,nres1)
c!!!	        stop ' gettotscore '
	      end if
              xtotscore=totscore
        end if
c
500	format(20i4)
c
        return
        end
c
c-------------------------------------------------------------------------------
c
	subroutine getblocks(ali1,ind,blocksize,block,nblock,nres1)
        implicit none
	include 'gagasizes.for'
	integer*2 ali1(maxres,popsize)
	integer blocksize,nblock,block(maxres),ind,nres1
c
	integer i,j

	do i=1,maxres
		block(i)=0
	end do
	nblock=0
	if(blocksize.gt.1) then
	  do i=1,nres1-blocksize+1
		if(ali1(i,ind).ne.0)then
			if(block(i).eq.0) nblock=nblock+1
			do j=0,blocksize-1
				block(i+j)=nblock
			end do
		end if
	  end do
	else
c	  blocksize=1 special case
	  if(ali1(1,ind).ne.0) then
		nblock=1
		block(1)=1
	  end if
	  do i=2,nres1
		if(ali1(i,ind).ne.0) then
c			special case: N-term of second sequence
			if(ali1(i,ind).eq.1) nblock=nblock+1
			if((ali1(i,ind).ne.ali1(i-1,ind)+1).and.
     $			   (ali1(i,ind).ne.ali1(i-1,ind)-1)) nblock=nblock+1
			block(i)=nblock
		end if
	  end do
	end if	

	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine readgenepool(filnam,ngene,genepool)
         implicit none
	include 'gagasizes.for'
	integer ngene
	character*80 filnam
	integer*2 genepool(5,maxpair)
	integer i,j

	open(90,file=filnam,status='old')
	read(90,500) ngene,((genepool(j,i),j=1,5),i=1,ngene)
	close(90)
500	format(10i8)
	write(*,*) ngene,' genes read in'

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine getali(ali1,a1,blocksize)
        implicit none
	include 'gagasizes.for'
	integer*2 ali1(maxres)
	integer blocksize,a1(maxres)
c
	integer j,l,i1

	do j=1,maxres
		a1(j)=0
	end do
	do j=1,maxres
		i1=ali1(j)
		if(i1.ne.0) then
			do l=0,blocksize-1
				a1(j+l)=i1+l
			end do
		end if
	end do
	
	return
	end
c
c	---------------------------------------------------------------------
c
	subroutine sortscores(d,score)
	implicit none
	include 'gagasizes.for'
	integer d(popsize),ve(popsize)
	real score(popsize)
	integer i
c
        do i=1,popsize
              d(i)=i
              ve(i)=nint(-score(i))
        end do
 	call j_index_Qsort(1,popsize,popsize,d,ve)
c
	return
	end
c
c	---------------------------------------------------------------------
c
