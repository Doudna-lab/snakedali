c This module/program is part of DaliLite (c) L. Holm 1999
c
	program gaga
c
c	gaga - main program
c
c	link gagarubber gagatool gagahash gagagene gagadccp clean dimple triplet
c		u3b
c
c	calculates a sequence-structure alignment under a pair potential
c	this program superimposes CA-CA distance maps
c
c	input: in setparameters & setup
c		d0 blocksize
c		dsspfile1
c		chainid1
c		dsspfile2
c		chainid2
c
c	parameters:
c		d0      -- distance-deviation threshold
c		blocksize --- fragment size
c
c		tolerance -- lfrozen score-change cutoff
c
c		ltop    -- impose same topology if .TRUE.
c		lpara	-- allow antiparallel matches if .TRUE.
c		lpretty -- pretty sequence alignments if .TRUE.
c
c		in gagasizes.for:
c		popsize -- how many alis are generated
c		maxd    -- max sum-of-distances in a square
c		pseed   -- only seed clusters with pairs scoring>pseed *maximal
c		maxpair -- how many pairs are saved
c		maxres  -- max size of protein
c		maxgen  -- max cycles in building up the alignment
c		maxhist -- max squares saved per sum-of-distances in a square
c
c	output:
c		code1code2.gene        -- fragment pair lists    (getgenes)
c		code1code2.ali         -- top popsize alignments (getalignments)
c		code2code2.dccp        -- neat output            (output)
c
	implicit none
	include 'gagasizes.for'
	character chainid1,chainid2,seq1(maxres0),seq2(maxres0),
     $		struc1(maxres0),struc2(maxres0)
	character*4 code1,code2,resno1(maxres0),resno2(maxres0)
	integer nres1,nres2,relacc1(maxres0),relacc2(maxres0)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	real ca1(3,maxres0),ca2(3,maxres0)
	character*80 hdr1,hdr2,infile
	character*80 dalidatpath_1, dalidatpath_2
	logical lfitz
c
	integer ierr
	logical file_exists
c
	integer*2 genepool(5,maxpair)
	integer ngene
	integer*2 ali1(maxres,popsize),preali1(maxres0)
	real score(popsize),oldscore,ds
	integer i,j,k,npr,pr1(2,1000),pr2(2,1000)
	logical ldef,lpreali,ltop0
	integer ntetra,nres0,tsil(maxres),nres02,tsil2(maxres)
	integer*2 tetrapool(2,maxpair)
c
	integer ali(maxres)
	real u(3,3),t(3),fitzrcut,rms
	integer nx,fitzmaxiter,nali,niter
	parameter (fitzrcut=4.0, fitzmaxiter=10)
c
	integer jnk1(dim1)
	common /junk1/jnk1
	integer jnk23(dim23)
	common /junk23/jnk23
	integer jnk4(dim4)
	common /junk4/jnk4
	integer jnk5(dim5)
	common /junk5/jnk5
	integer jnk6(dim6)
	common /junk6/jnk6
c
	ldef=.true.
	open(90,file='dali.default',status='old',err=5)
	goto 6
5	ldef=.false.
	write(*,*) ' *** dali.default not found *** that''s ok ***'
	goto 7
6	call defaults(90,ldef)
	write(*,*),'defaults done'
	close(90)
7	continue
	
	write(*,*) 'enter dalidatpath_1'
	read(*,501) dalidatpath_1
	write(*,*) 'enter dalidatpath_2'
	read(*,501) dalidatpath_2
501	format(a80)
	write(*,*) 'enter lfitz'
	read(*,*) lfitz
c
	write(*,*) ' fitz switch is: ',lfitz
	write(*,*) ' topology switch is: ',ltop
	if((.not.ltop).and.lpara) write(*,520) 'allowed'
	if(ltop.or.(.not.lpara)) write(*,520) 'disallowed'
	call weights
1	do i=1,maxres
		tsil(i)=0
	end do
	call setup(code1,chainid1,nres1,seq1,struc1,relacc1,resno1,ca1,d1,
     $		hdr1,ierr,infile,lpreali,pr1,pr2,npr,nres0,tsil,
     $		dalidatpath_1)
	if(ierr.eq.-99) stop ' END -- Normal termination '
	if(ierr.eq.-1) then
		write(*,*) ' file not found -- return to top level '
		goto 1
	end if

	write(*,*) code1,' query ', nres1

c
c	main loop: exit on file-not-found
c
10	do i=1,maxres
		tsil2(i)=0
	end do
	call setup(code2,chainid2,nres2,seq2,struc2,relacc2,resno2,ca2,d2,
     $		hdr2,ierr,infile,lpreali,pr1,pr2,npr,nres02,tsil2,
     $		dalidatpath_2)
	if(ierr.eq.-99) then
		write(*,*) ' END of comparison list -- return to top level '
c!		close(91)
		goto 1
	end if
	if(ierr.eq.-1) then
		write(*,*) ' file not found -- comparison skipped '
		goto 10
	end if
	write(*,*) code2,' target ',nres2
	write(*,500) code1,chainid1,nres1,code2,chainid2,nres2,lpreali
c
c	convert prealignment ranges pr1,pr2 to preali1
c
	ltop0=ltop
	call initali1(preali1,nres1)
	write(*,*),'initali done',nres1,nres2,nres0
	nx=0
	if(lpreali) then
		if(.not.ltop) then
			write(*,*) 'prealignment implies sequential constraint'
			ltop=.true.
		end if
		do i=1,npr
			k=pr2(1,i)
			do j=pr1(1,i),pr1(2,i)
				if(j.le.maxres.and.k.le.maxres) preali1(j)=k
				k=k+1
				nx=nx+1
			end do
		end do
29	end if
	write(*,*),'preali1 done',k,(tsil(j),j=1,nres1),
     $          (tsil2(j),j=1,nres2)
c
c	run fitz to extend short prealignment
c
	if(lfitz) then
!		write(*,*) 'initial preali1',nres1,(preali1(i),i=1,nres1)
c		preali1 is integer*2, getut expects ali as integer
		do i=1,nres1
			ali(i)=preali1(i)
		end do
c		initial superimposition by u3b
		call getut(nres1,ali,ca1,ca2,u,t,nali,rms)
!		write(*,*) 'getut result:',nali,rms,nres1,(ali(i),i=1,nres1)
		call transrotate(ca1,nres1,u,t)
c		fitz iterations
		call fitz(ca1,nres1,ca2,nres2,ali,fitzrcut,fitzmaxiter,rms,
     $                  nali,niter)
c		overwrite preali1 if fitzed is longer
		if(nali.gt.nx) then 
			do i=1,nres1
				preali1(i)=ali(i)
!			write(*,*) 'fitz result ',i,ali(i),preali1(i)
			end do
			write(*,*) 'fitzed preali1',rms,nali,niter,
     $                          (preali1(i),i=1,nres1)
		end if
	end if
c
c	reverse Polish stack mode of operation
c
c	if(file_exists(pairfilnam(code1,chainid1,code2,chainid2,'gene'))) then
c	    call readgenepool(pairfilnam(code1,chainid1,code2,chainid2,
c     $		'gene'),ngene,genepool)
c	else if(lpreali) then
	    call compresspreali(preali1,tsil,tsil2,nres1,nres2,nres0)
	    write(*,*),'compresspreali done'
	    call testi(nres1,nres2,preali1,d1,d2,ntetra,tetrapool)
	    write(*,*),'testi done',nres1,nres2,ntetra
	    lverb=.false.
	    score(1)=0.0
	    do i=1,maxres
		ali1(i,1)=0
	    end do
	    if(ntetra.gt.0) call lean_mc(50,4,nres1,nres2,d1,d2,.true.,
     $		preali1,itrim(2),score,ali1,ntetra,tetrapool)
	    write(*,*),'lean_mc done'
c
c	    iterate until improvement < 10
c
	    oldscore=score(1)
	    ds=oldscore
	    do while(ds.gt.refitol)
	      do i=1,maxres
		preali1(i)=ali1(i,1)
	      end do
	      call testi(nres1,nres2,preali1,d1,d2,ntetra,tetrapool)
	      lverb=.false.
	      score(1)=0.0
	      do i=1,maxres
		ali1(i,1)=0
	      end do
	      if(ntetra.gt.0) call lean_mc(50,4,nres1,nres2,d1,d2,.true.,
     $		preali1,itrim(2),score,ali1,ntetra,tetrapool)
	      ds=score(1)-oldscore
	      oldscore=score(1)
	    end do
	    write(*,*) (ali1(k,1),k=1,nres1)
c
	    write(*,*) 'call output'
	    call output(code1,chainid1,code2,chainid2,nres1,nres2,seq1,
     $	seq1,struc1,struc2,resno1,resno2,1,hdr1,hdr2,ca1,ca2,
     $	relacc1,relacc2,d1,d2,ali1,score,nres0,tsil,nres02,tsil2)
	    goto 210
c	else
c	    call getgenes(code1,chainid1,code2,chainid2,d1,d2,nres1,nres2,
c     $		geneblocksize,ngene,genepool,lpreali,preali1) !!! bugged !!!
c	end if
	do j=1,nrepeat
	    do i=1,popsize
		score(i)=0.0
	    end do
!	    write(*,*),'getali: ',code1,chainid1,code2,chainid2,nres1,nres2,ngene
	    call getalignments(code1,chainid1,code2,chainid2,
     $		d1,d2,nres1,nres2,geneblocksize,outblocksize,
     $		ngene,genepool,ali1,score,lpreali,preali1)
c
!	    write(*,*),'output: ',code1,chainid1,code2,chainid2,nres1,nres2
	    call output(code1,chainid1,code2,chainid2,nres1,nres2,seq1,
     $	seq2,struc1,struc2,resno1,resno2,outblocksize,hdr1,hdr2,ca1,ca2,
     $	relacc1,relacc2,d1,d2,ali1,score,nres0,tsil,nres02,tsil2)
	end do
210	ltop=ltop0
	goto 10
c
c	end of main loop
c
500	format(/' comparing: ',a4,a1,i5,' res to ',a4,a1,i5,' res ',l5)
510	format(a20,' output file exists already !!!')
520	format(' antiparallel matches are ',a20/)

	end
c
c----------------------------------------------------------------------
c
	subroutine compresspreali(preali1,tsil,tsil2,nres1,nres2,nres0)
	implicit none
	include 'gagasizes.for'
	integer*2 preali1(maxres0)
	integer tsil(maxres),nres1,tsil2(maxres),nres2,nres0
c
	integer list(maxres0),list2(maxres0),i,k,tmp(maxres0)
c
	write(*,*),'this is compresspreali ',nres1,nres2,nres0,maxres,
     $          maxres0
	do i=1,maxres0
		list(i)=0
		list2(i)=0
		tmp(i)=0
	end do
	do i=1,min(maxres0,nres1)
		list(tsil(i))=i
	end do
	do i=1,min(maxres0,nres2)
		list2(tsil2(i))=i
	end do
	do i=1,min(maxres0,nres0)
		k=preali1(i)
		if(k.gt.0) then
			if(list(i).gt.0) tmp(list(i))=list2(k)
		end if
	end do
	do i=1,maxres0
		preali1(i)=tmp(i)
	end do
c
	write(*,*) 'compressed preali1:',(preali1(i),i=1,nres1)
c
	return
	end
c
c----------------------------------------------------------------------
c
	function scorefun(a,b)
	implicit none
	include 'gagasizes.for'
	real scorefun
	integer*2 a,b
c
	real x,y
        integer yint
!       hardcoded parameters enveloperadius=20.0 and d0=0.20 - LH 27 Jun 05
c !!! 	elastic uses weights !!!
	x=float(abs(a-b))/10.0
!	if(lela) then
		y=float(a+b)
                y=y/20.0
                yint=nint(y)
		if(y.gt.100) then
			scorefun=0.0
		else
			if(y.gt.0) then
				scorefun=weight(yint)*(0.20-x/y)
			else
				scorefun=weight(yint)*0.20
			end if
		end if
		scorefun=scorefun*100
!	else
c		rigid
!		scorefun=d0-x
!		scorefun=scorefun*100
!	end if
        !
!        write(*,*) 'scorefun',a,b,x,y,yint,x/y,scorefun
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine weights
	implicit none
	include 'gagasizes.for'
	integer i
c
	real x
c
	x=1/(enveloperadius*enveloperadius) ! 1/20.0/20.0=0.0025
	do i=0,100
                x=-0.0025*i*i
		weight(i)=exp(x) 
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine defaults(iunit,ldef)
	implicit none
	include 'gagasizes.for'
c
	integer nparam,iunit
	logical ldef
c
	parameter(nparam=48)
	character*20 key(nparam)
	data key/'machine','dummy','dalidatpath',
     $		'd0','enveloperadius','tripletscorecutoff',
     $		'killcutoff','killcutoff','killcutoff',
     $		'pblock','corecutoff','fokcutoff','putcutoff',
     $		'plimit','p600','p720',
     $		'iounit1','geneblocksize','outblocksize','maxstep',
     $		'itrim','itrim','itrim',
     $		'ikill','ikill','ikill','nclone','nclone','nclone',
     $		'ltop','lpara','lali','lpretty','lverb','seed',
     $		'nrepeat','lsavegenes','nbest','printpopsize','lcalpha',
     $		'seedcutoff','refitol','hack0','hack1','lplottriplets',
     $		'ldebug','lprintparams','ltube'/
c
	character*20 word
	character c
	integer i,j
	character*80 dalidatpath
c
c	set defaults
c
	machine='U'
	dalidatpath='/junk/holm/tmp1/'
	d0=0.2
	enveloperadius=20.0
	tripletscorecutoff=0.0
	killcutoff(1)=0.50
	killcutoff(2)=0.80
	killcutoff(3)=1.00
	pblock=0.70
	corecutoff=20
	fokcutoff=-1000.0
	putcutoff=0.5
	plimit=900
	p600=900
	p720=1000
	iounit1=90
	geneblocksize=6
	outblocksize=4
	maxstep=100
	itrim(1)=1
	itrim(2)=5
	itrim(3)=5
	ikill(1)=1
	ikill(2)=20
	ikill(3)=20
	nclone(1)=1
	nclone(2)=1
	nclone(3)=popsize
	ltop=.false.
	lpara=.true.
	lali=.false.
	lpretty=.false.
	lverb=.false.
	seed=124462287
	nrepeat=1
	lsavegenes=.false.
	nbest=1
	printpopsize=nbest
	lcalpha=.false.
	seedcutoff=0.10
	refitol=100.0
	hack0=100.0
	hack1=1.0
	lplottriplets=.false.
	ldebug=.false.
	lprintparams=.false.
	ltube=.false.
c
c
c
	if(ldef) then
10	  read(iunit,500,end=19) c,word
		if(c.eq.'#') then
			i=0
			do while(i.lt.nparam)
				i=i+1
				if(word.eq.key(i)) then
		if(i.eq.1) read(iunit,*) machine
		if(i.eq.3) read(iunit,510) dalidatpath
		if(i.eq.4) read(iunit,*) d0
		if(i.eq.5) read(iunit,*) enveloperadius
		if(i.eq.6) read(iunit,*) tripletscorecutoff
		if(i.eq.7) read(iunit,*) (killcutoff(j),j=1,3)
		if(i.eq.10) read(iunit,*) pblock
		if(i.eq.11) read(iunit,*) corecutoff
		if(i.eq.12) read(iunit,*) fokcutoff
		if(i.eq.13) read(iunit,*) putcutoff
		if(i.eq.14) read(iunit,*) plimit
		if(i.eq.15) read(iunit,*) p600
		if(i.eq.16) read(iunit,*) p720
		if(i.eq.17) read(iunit,*) iounit1
		if(i.eq.18) read(iunit,*) geneblocksize
		if(i.eq.19) read(iunit,*) outblocksize
		if(i.eq.20) read(iunit,*) maxstep
		if(i.eq.21) read(iunit,*) (itrim(j),j=1,3)
		if(i.eq.24) read(iunit,*) (ikill(j),j=1,3)
		if(i.eq.27) read(iunit,*) (nclone(j),j=1,3)
		if(i.eq.30) read(iunit,*) ltop
		if(i.eq.31) read(iunit,*) lpara
		if(i.eq.32) read(iunit,*) lali
		if(i.eq.33) read(iunit,*) lpretty
		if(i.eq.34) read(iunit,*) lverb
		if(i.eq.35) read(iunit,*) seed
		if(i.eq.36) read(iunit,*) nrepeat
		if(i.eq.37) read(iunit,*) lsavegenes
		if(i.eq.38) read(iunit,*) nbest
		if(i.eq.39) read(iunit,*) printpopsize
		if(i.eq.40) read(iunit,*) lcalpha
		if(i.eq.41) read(iunit,*) seedcutoff
		if(i.eq.42) read(iunit,*) refitol
		if(i.eq.43) read(iunit,*) hack0
		if(i.eq.44) read(iunit,*) hack1
		if(i.eq.45) read(iunit,*) lplottriplets
		if(i.eq.46) read(iunit,*) ldebug
		if(i.eq.47) read(iunit,*) lprintparams
		if(i.eq.48) read(iunit,*) ltube
				i=nparam
				end if
			end do
		end if
	  goto 10
19	end if
c
c	check parameters make sense
c
	if((machine.ne.'U').and.(machine.ne.'V')) stop 'unknown machine'
	if(outblocksize.gt.geneblocksize) stop 'illegal blocksizes'
	if(nbest.gt.popsize) then
		write(*,*) ' resetting nbest to popsize'
		nbest=popsize
	end if
	if(printpopsize.gt.popsize) then
		write(*,*) ' resetting printpopsize to popsize'
		printpopsize=popsize
	end if
	if(ltube.and.(.not.ltop)) then
		write(*,*) ' ltube option on, so using sequential constraint '
		ltop=.false.
	end if
c
c	automatically select elastic/rigid score
c

	lela=(d0.lt.1.0)
	if(lela) write(*,*) ' using elastic score'
	if(.not.lela) write(*,*) ' using rigid score'
	write(*,*) ' refine scores > ',hack0,' or ',hack1,'*nres1'
c
c	echo parameter settings
c
	if(lprintparams) then
	  do i=1,nparam
		if(i.eq.1) write(*,*) key(i),': ', machine
		if(i.eq.3) write(*,*) key(i),': ', dalidatpath
		if(i.eq.4) write(*,*) key(i),': ', d0
		if(i.eq.5) write(*,*) key(i),': ', enveloperadius
		if(i.eq.6) write(*,*) key(i),': ', tripletscorecutoff
		if(i.eq.7) write(*,*) key(i),': ', (killcutoff(j),j=1,3)
		if(i.eq.10) write(*,*) key(i),': ', pblock
		if(i.eq.11) write(*,*) key(i),': ', corecutoff
		if(i.eq.12) write(*,*) key(i),': ', fokcutoff
		if(i.eq.13) write(*,*) key(i),': ', putcutoff
		if(i.eq.14) write(*,*) key(i),': ', plimit
		if(i.eq.15) write(*,*) key(i),': ', p600
		if(i.eq.16) write(*,*) key(i),': ', p720
		if(i.eq.17) write(*,*) key(i),': ', iounit1
		if(i.eq.18) write(*,*) key(i),': ', geneblocksize
		if(i.eq.19) write(*,*) key(i),': ', outblocksize
		if(i.eq.20) write(*,*) key(i),': ', maxstep
		if(i.eq.21) write(*,*) key(i),': ', (itrim(j),j=1,3)
		if(i.eq.24) write(*,*) key(i),': ', (ikill(j),j=1,3)
		if(i.eq.27) write(*,*) key(i),': ', (nclone(j),j=1,3)
		if(i.eq.30) write(*,*) key(i),': ', ltop
		if(i.eq.31) write(*,*) key(i),': ', lpara
		if(i.eq.32) write(*,*) key(i),': ', lali
		if(i.eq.33) write(*,*) key(i),': ', lpretty
		if(i.eq.34) write(*,*) key(i),': ', lverb
		if(i.eq.35) write(*,*) key(i),': ', seed
		if(i.eq.36) write(*,*) key(i),': ', nrepeat
		if(i.eq.37) write(*,*) key(i),': ', lsavegenes
		if(i.eq.38) write(*,*) key(i),': ', nbest
		if(i.eq.39) write(*,*) key(i),': ', printpopsize
		if(i.eq.40) write(*,*) key(i),': ', lcalpha
		if(i.eq.41) write(*,*) key(i),': ', seedcutoff
		if(i.eq.42) write(*,*) key(i),': ', refitol
		if(i.eq.43) write(*,*) key(i),': ', hack0
		if(i.eq.44) write(*,*) key(i),': ', hack1
		if(i.eq.45) write(*,*) key(i),': ', lplottriplets
		if(i.eq.46) write(*,*) key(i),': ', ldebug
		if(i.eq.47) write(*,*) key(i),': ', lprintparams
		if(i.eq.48) write(*,*) key(i),': ', ltube
	  end do
	end if
c
500	format(a1,a20)
510	format(a80)
c
	return
	end

