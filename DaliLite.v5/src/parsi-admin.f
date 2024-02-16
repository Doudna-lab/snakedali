c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
c	goal: align on any pairwise property
c	- Daliscore
c	- contact strength
c	- other energy term
c	- solvent accessibility
c	- fosfos contact preference
c	- other contact preference (Sippl, ...)
c	- other singleton terms (PHD, ...)
c	+ weighted combinations
c	+ weighted profiles
c
c	affected subroutines:
c	- getstructure1
c	- getsequence2
c	- segsegscore
c	- setcut
c	- weights
c	- fillscoretable
c
c----------------------------------------------------------------------
c
c	administration module
c
c----------------------------------------------------------------------
c
	subroutine setcut(cut,ndom,domns,domseglist,segmentrange)
	implicit none
	include 'parsizes.for'
	integer cut(maxdom),ndom,domns(maxdom),domseglist(maxseg,maxdom)
	integer segmentrange(2,maxseg)
c
	integer idom,lali,seglist(maxdom),getlali,i
	real x,l
c
	do idom=1,ndom
		do i=1,domns(idom)
			seglist(i)=domseglist(i,idom)
		end do
		lali=getlali(segmentrange,domns(idom),seglist)
c!		cut(idom)=210*lali*lali
c		lali vs. score
		l=min(float(lali),200.0)
c
c		mean minus one sigma
c		x=10000.0*(0.83259+0.11186*lali+3.3537e-5*lali*lali*lali
c     $			+1.475e-3*lali*lali-1.579e-7*lali*lali*lali*lali
c     $			-1.2956+1.6648e-2*lali-9.945e-4*lali*lali)
c
c		90 % of mean
		x=9000.0*(0.83259+0.11186*l+3.3537e-5*l*l*l
     $			+1.475e-3*l*l-1.579e-7*l*l*l*l)
		cut(idom)=max(x,0.0)
c		cut(idom)=0
c		write(*,*) idom,cut(idom)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine setngap(ngap,segmentrange,minseglen,bl,nseg)
	implicit none
	include 'parsizes.for'
	integer ngap(maxseg),minseglen(maxseg),segmentrange(2,maxseg)
	integer nseg,bl
c
	integer lres,iseg,a1,a2,minlen
c
	do iseg=1,nseg
		a2=segmentrange(2,iseg)
		a1=segmentrange(1,iseg)
		minlen=minseglen(iseg)

		lres=(a2-a1+1-minlen)
		if(lres.lt.-29) then
			write(*,*) ' forcing long gap from',lres,' to -29'
			lres=-29
		end if
		ngap(iseg)=(lres+bl-1)/bl

	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine setldom(ldom,ndom,node_size,startsize)
c
c	ldom(idom)=.true. means align this node
c
	implicit none
	include 'parsizes.for'
	logical ldom(0:maxdom)
	integer ndom,node_size(maxdom),startsize
c
	integer idom
c
	ldom(0)=.false.
	do idom=1,ndom
		ldom(idom)=(node_size(idom).ge.startsize)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine nextprotein(i,string,text,ldb,curprot,lastprot,
     $ list,flag)
c
c	flag=0: ok
c	flag=1: end
c
	implicit none
	include 'parsizes.for'
	character*10 string
	integer i,flag
	character*5 list(maxprot),cd
	character*10 text
	logical ldb
	integer curprot,lastprot
c
	flag=0
	if(ldb) then
		curprot=curprot+1
		if(curprot.gt.lastprot) flag=1
		if(flag.eq.0) cd=list(curprot)
	else
  		write(*,*) ' enter ',text,' code+chainid (END to quit) ?'
 		read(*,500) cd
  		if(cd(1:3).eq.'END') flag=1
	end if
	if(cd(5:5).eq.' ') cd(5:5)='_'
	if(flag.eq.0) string(i:i+4)=cd
c
500     format(a5)
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getstructure1(nseg,segmentrange,na1,nb1,na2,nb2,
     $		ndom,node_child,dist,nres,string,upper,lower,ss,checkrange,
     $		secstr,lfix,lfix1,cut,domns,domseglist,flag,checkx,dist1sum,
     $		segmentrange0,ldb,dalidatpath)
	implicit none
	include 'parsizes.for'
	character*80 dalidatpath
	integer nres,nseg,segmentrange(2,maxseg),na1,nb1,ndom,na2,nb2
	integer node_child(2,maxdom),upper(maxseg,maxseg),
     $ lower(maxseg,maxseg)
	integer*2 dist(maxres1*maxres1)
	character*10 string
	integer ss(maxseg,maxseg),flag,cut(maxdom),segmentrange0(2,maxseg)
	logical lfix(maxseg,maxdom)
	integer domns(maxdom)
	integer domseglist(maxseg,maxdom)
	logical lfix1(maxseg,maxdom),ldb
	integer checkrange(2,maxseg),checkx(maxseg),
     $ dist1sum(maxres1,0:maxres1)
	character secstr(maxseg),node_type(maxdom)
c
	real ca(3,maxres)
	character chainid
	character*4 code
c
	nseg=0
	na1=0
	nb1=0
	ndom=0
	nres=0
	flag=0
	code=string(1:4)
	chainid=string(5:5)
c	write(*,*) 'getprotein1 ',code,chainid,' ',dalidatpath
	call getprotein1(code,chainid,90,ca,nres,
     $		segmentrange0,nseg,na1,nb1,ndom,node_child,checkrange,
     $		secstr,checkx,node_type,domns,domseglist,ldb,dalidatpath)
c	write(*,*) 'getprotein1 done',nseg,ndom,nseg,na1,na2,nb1,nb2
	if((nseg.le.2).or.(ndom.lt.nseg).or.(na1*na2+nb1*nb2.eq.0))then
                flag=1
                return
        end if 
c	write(*,*) 'calling compressca'
	call compressca(ca,nres,nseg,segmentrange0,segmentrange)
c	write(*,*) 'compressca done'
	if(lskipflag) then
                flag=2
                return
        end if
c	write(*,*) 'calling getdist'
	call getdist(nres,ca,dist)
c	write(*,*) 'getdist done'
	call selfscore(nseg,segmentrange,dist,nres,ss)
	call getupperlower(dist,nres,nseg,segmentrange,lower,upper)
c
c	hack: set upper distance limit to 400
c
	call hackdist(dist,nres)
c
	call flex(ss,nseg,ndom,domns,domseglist,lfix,lfix1)
c	write(*,*) 'flex done'
	call setcut(cut,ndom,domns,domseglist,segmentrange)
	call getdist2sum(nres,ca,dist1sum)
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine setminseglen(secstr,minseglen,nseg,lene,lenh)
	implicit none
	include 'parsizes.for'
	character secstr(maxseg)
	integer minseglen(maxseg),nseg,lene,lenh
c
	integer iseg
c
	do iseg=1,nseg
		minseglen(iseg)=lene
		if(secstr(iseg).eq.'H') minseglen(iseg)=lenh
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getsequence2(nres,na2,nb2,dist2,string,dist2sum,
     $		segment,nseg,secstr,na1,nb1,flag,ldb,dalidatpath)
	implicit none
	include 'parsizes.for'
	character*80 dalidatpath
	integer nres,na2,nb2
	integer*2 dist2(maxres2*maxres2)
	character*10 string
	integer dist2sum(maxres2,0:maxres2),flag,na1,nb1
	integer segment(maxres2),nseg
	character secstr(maxseg)
	logical ldb
c
	real ca(3,maxres)
c
	nres=0
	na2=0
	nb2=0
	flag=0
        lskipflag=.false. ! may be set to .true. later
	call getprotein2(string(6:9),string(10:10),90,ca,nres,
     $		na2,nb2,segment,nseg,secstr,ldb,dalidatpath)
	if((nres.eq.0).or.(na1*na2+nb1*nb2.eq.0)) flag=1
c	write(*,*) 'flag=',flag
	if(flag.ne.0) return
	if(nres.gt.maxres2) then
	write(*,*)  'FATAL ERROR -- nres overflow: getsequence'
		lskipflag=.true.
                flag=2
		return
	end if

	call getdist2sum(nres,ca,dist2sum)
c
c	hack distance matrix to bin-indices
c
	call getdist2(nres,ca,dist2)
c
	return
	end
c
c----------------------------------------------------------------------
c
        subroutine weights
        implicit none
        include 'parsizes.for'
        integer i
c
        real x,enveloperadius
	parameter(enveloperadius=20.0)
c
        x=1/(enveloperadius*enveloperadius)
        do i=1,1000
                weight(i)=nint(100*exp(-x/100*i*i))
		if(weight(i).lt.5) weight(i)=0
c		if(mod(i,40).eq.0) write(*,500) i,weight(i)
        end do

500	format('weight:',i5,i10)

        return
        end
c
c----------------------------------------------------------------------
c
	subroutine fillscoretable
	implicit none
	include 'parsizes.for'
c
	integer a,b
	real x,y
c
	do b=1,160
		if(b.lt.100) then
			y=b*0.1
		else if(b.lt.125) then
			y=10.0+(b-100)*0.4
		else
			y=b-125+20.0
		end if
		do a=1,400
			x=a/10.0
			scoretable(b,a)=nint(100.0*weight(a)*(0.20-abs(x-y)/x))
c	if(scoretable(b,a).gt.0) write(95,*) a,b,x,y,scoretable(b,a)
		end do
	end do

        return
        end
c
c------------------------------------------------------------------------------
c
	function getlali(segmentrange,ns,seglist)
	implicit none
	include 'parsizes.for'
	integer getlali,segmentrange(2,maxseg),ns,seglist(maxseg)
c
	integer j,lali
c
	lali=0
	do j=1,ns
		lali=lali+segmentrange(2,seglist(j))
     $			 -segmentrange(1,seglist(j))+1
	end do
c
c	set upper limit on domain size = 200 residues !
c
	getlali=min(lali,200)

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getprotein1(code,chainid,iunit,ca,nres,segmentrange,
     $ nseg,na,nb,ndom,node_child,checkrange,secstr,checkx,
     $ node_type,domns,domseglist,ldb,dalidatpath)
	implicit none
	include 'parsizes.for'
	character*80 dalidatpath
	character*4 code
	character chainid
	integer nres,nseg,iunit
	character secstr(maxseg)
	integer segmentrange(2,maxseg),checkrange(2,maxseg),checkx(maxseg)
	real ca(3,maxres)
	integer ndom,node_child(2,maxdom),na,nb
	character node_type(maxdom)
	integer domns(maxdom),domseglist(maxseg,maxdom)
	logical ldb
c	logical ldb was to prevent type mismatch error(from integer) JONG
c
	character*80 filnam,constructfilnam
	character*5 cd
c
	cd(1:4)=code
	cd(5:5)=chainid
	filnam=constructfilnam(cd,dalidatpath,'.dat')
c	write(*,*) 'read file ',filnam
10	open(iunit,file=filnam,status='old',err=99)
	call parsireadproteindata(nres,nseg,na,nb,secstr,segmentrange,
     $		checkrange,checkx,ca,ndom,node_type,node_child,domns,
     $		domseglist,iunit)
	close(iunit)
	call treehack(ndom,node_child,domns,domseglist)
c	write(*,*) 'treehack done'
c
	return
c	error exit
99	write(*,*) 'file not found ',filnam
!	if(ldb) then
!		write(*,*) 'type CONTINUE to try again'
!		goto 10
!	end if
	nres=0

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine compressca(ca,nres1,nseg,segmentrange0,segmentrange)
	implicit none
	include 'parsizes.for'
	integer nseg,segmentrange(2,maxseg),segmentrange0(2,maxseg)
	integer nres1
	real ca(3,maxres)
c
	integer l,iseg,sr(2,maxseg),ires,a,i
	real cx(3,maxres)
c
c	shortlist
c
	l=0
        lskipflag=.false. ! may be set to .true. later
	do iseg=1,nseg
		l=l+1
		sr(1,iseg)=l
		l=l+segmentrange0(2,iseg)-segmentrange0(1,iseg)
		sr(2,iseg)=l
	end do
	do iseg=1,nseg
		do ires=sr(1,iseg),sr(2,iseg)
			if(ires.gt.maxres1) then
	write(*,*)  'FATAL ERROR -- maxres1 overflow: compressca ',ires
				lskipflag=.true.
				return
			end if
			a=segmentrange0(1,iseg)+ires-sr(1,iseg)
			do i=1,3
				cx(i,ires)=ca(i,a)
			end do
		end do
	end do
c	write(*,*) 'length original:',nres1,' compressed:',l
c
c	overwrite nres1,segmentrange,dist
c	keep original segmentrange !
c
	nres1=sr(2,nseg)
	do ires=1,nres1
		do i=1,3
			ca(i,ires)=cx(i,ires)
		end do
	end do
	do iseg=1,nseg
		segmentrange(1,iseg)=sr(1,iseg)
		segmentrange(2,iseg)=sr(2,iseg)
	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
        subroutine getdist(nres1,ca,dist)
        implicit none
        include 'parsizes.for'
c
c	fills dist(nres,nres) !!
c
        integer nres1
        real ca(3,nres1)
	integer*2 dist(nres1,nres1)
c
        integer*2 parsidistance
c
        integer i,j
c
        do i=1,nres1
c		!minimal distance=1, maximal distance=maxdist
                dist(i,i)=1
                do j=i+1,nres1
                  dist(i,j)=parsidistance(ca(1,i),ca(2,i),
     $                  ca(3,i),ca(1,j),ca(2,j),ca(3,j))
		  if(dist(i,j).gt.1000) dist(i,j)=1000
                  dist(j,i)=dist(i,j)
                end do
        end do

        return
        end
c
c------------------------------------------------------------------------------
c
	subroutine selfscore(nseg,segmentrange,dist,nres1,ss)
        implicit none
        include 'parsizes.for'
	integer nres1
        integer*2 dist(nres1,nres1)
        integer nseg,segmentrange(2,maxseg)
        integer ss(maxseg,maxseg)
c
        integer iseg,jseg,ires,jres
        integer s,p,stot,lali
c
c       calculate segment-segment self-Dali-scores
c               use as maximum possible seg-seg estimate
c
        do iseg=1,nseg
	  do jseg=1,iseg
          s=0
          do ires=segmentrange(1,iseg),segmentrange(2,iseg)
           do jres=segmentrange(1,jseg),segmentrange(2,jseg)
		p=dist(ires,jres)
		s=s+weight(p)*20
           end do
          end do
c
c	  cutoff 10 % of db-mean, i.e. 10 % of 420*lali*lali
c
	  lali=segmentrange(2,iseg)-segmentrange(1,iseg)+1+
     $		segmentrange(2,jseg)-segmentrange(1,jseg)+1
	  if(s.lt.42*lali*lali) s=0
          ss(iseg,jseg)=s
	  ss(jseg,iseg)=s
         end do
        end do
	stot=0
	do iseg=1,nseg
	 do jseg=1,nseg
		stot=stot+ss(iseg,jseg)
	 end do
c	 write(*,*) 'seg:',iseg,(ss(iseg,jseg),jseg=1,nseg)
	end do
c	write(*,*) 'selfscore:',stot

        return
        end
c
c----------------------------------------------------------------------
c
	subroutine getupperlower(dist,nres1,nseg,segmentrange,lower,upper)
	implicit none
	include 'parsizes.for'
	integer nseg,segmentrange(2,maxseg)
	integer nres1,lower(maxseg,maxseg),upper(maxseg,maxseg)
	integer*2 dist(nres1,nres1)
c
	integer ires,jres,s,iseg,jseg
c
c	calculate 60 % and 140 % distance-sum bounds per segment-segment blocks
c
	do iseg=1,nseg
		do jseg=1,iseg
			s=0
			do ires=segmentrange(1,iseg),segmentrange(2,iseg)
			  do jres=segmentrange(1,jseg),segmentrange(2,jseg)
				s=s+dist(ires,jres)
			  end do
			end do
			lower(iseg,jseg)=70*s
			lower(jseg,iseg)=lower(iseg,jseg)
			upper(iseg,jseg)=130*s
			upper(jseg,iseg)=upper(iseg,jseg)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine hackdist(dist,nres1)
	implicit none
	include 'parsizes.for'
	integer nres1
	integer*2 dist(nres1,nres1)
c
	integer i,j
c
	do i=1,nres1
		do j=1,i-1
			if(dist(i,j).gt.400) then
				dist(i,j)=400
				dist(j,i)=400
			end if
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine flex(ss,nseg,ndom,domns,domseglist,lfix,lfix1)
	implicit none
	include 'parsizes.for'
	integer ss(maxseg,maxseg),nseg
	integer ndom,domns(maxdom),domseglist(maxseg,maxdom)
	logical lfix(maxseg,maxdom),lfix1(maxseg,maxdom)
c
	integer iseg,jseg,si(maxseg),is,idom,ns,seglist(maxseg),n,i
	real sf(maxseg,maxseg)
c
	do iseg=1,nseg
		si(iseg)=0
		do jseg=1,nseg
			if(jseg.ne.iseg) si(iseg)=si(iseg)+ss(iseg,jseg)
			sf(iseg,jseg)=0.0
		end do
		if(si(iseg).eq.0) si(iseg)=1
		do jseg=1,nseg
		  if(jseg.ne.iseg) sf(iseg,jseg)=float(ss(iseg,jseg))/si(iseg)
		end do
	end do
c
c	print free/fixed segments per node
c
	do idom=ndom,1,-1
		ns=domns(idom)
		do i=1,ns
			seglist(i)=domseglist(i,idom)
		end do
		do iseg=1,nseg
			lfix(iseg,idom)=.false.
		end do
		call getlfix(idom,lfix,0.70,ns,seglist,sf)
		call getlfix(idom,lfix1,0.90,ns,seglist,sf)
c
c		if domain-size.gt.2 fix densest segment if none fixed by cutoff !
c
		if(ns.gt.2) then
			n=0
			do is=1,ns
c				lfix1(seglist(is),idom)=lfix(seglist(is),idom)
				if(lfix(seglist(is),idom)) n=n+1
			end do
			if(n.eq.0) then
c!
				lfix(seglist(1),idom)=.true.
			end if
		end if
	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getprotein2(code,chainid,iunit,ca,nres,na,nb,segment,
     $		nseg,secstr,ldb,dalidatpath)
	implicit none
	include 'parsizes.for'
	character*4 code
	character chainid
	integer nres,nseg,iunit
	character secstr(maxseg)
	integer segmentrange(2,maxseg),segment(maxres2)
	real ca(3,maxres)
	integer ndom,node_child(2,maxdom),na,nb
	logical ldb
	character*80 dalidatpath,constructfilnam
c
	character*80 filnam
	integer ires,iseg,domns(maxdom),domseglist(maxseg,maxdom)
	integer checkrange(2,maxseg),checkx(maxseg)
	character node_type(maxdom)
	character*5 cd
c
	cd(1:4)=code
	cd(5:5)=chainid
	filnam=constructfilnam(cd,dalidatpath,'.dat')
c	write(*,*) 'read file ',filnam
	open(iunit,file=filnam,status='old',err=19)
	goto 20
c	reserve: look in CWD
	filnam=constructfilnam(cd,'./','.dat')
19	open(iunit,file=filnam,status='old',err=99)
20	call parsireadproteindata(nres,nseg,na,nb,secstr,segmentrange,
     $		checkrange,checkx,ca,ndom,node_type,node_child,domns,
     $		domseglist,iunit)
	close(iunit)
c
	do ires=1,nres
		segment(ires)=0
	end do
	do iseg=1,nseg
		do ires=segmentrange(1,iseg),segmentrange(2,iseg)
			segment(ires)=iseg
		end do
	end do
c
	return
c	error exit
99	write(*,*) 'file not found ',filnam
!	if(ldb) then
!		write(*,*) 'type CONTINUE to try again'
!		pause
!	end if
	nres=0

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getdist2sum(nres2,ca,dist2sum)
	implicit none
	include 'parsizes.for'
	integer nres2,dist2sum(nres2,0:nres2)
	real ca(3,maxres)
c
	integer ires,jres,distanceint4
c
	do ires=1,nres2
		dist2sum(ires,0)=0
		do jres=1,nres2
		  dist2sum(ires,jres)=dist2sum(ires,jres-1)+
     $			distanceint4(ca(1,ires),ca(2,ires),ca(3,ires),
     $			ca(1,jres),ca(2,jres),ca(3,jres))
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getdist2(nres2,ca,dist2)
        implicit none
        include 'parsizes.for'
        integer nres2
        real ca(3,nres2)
	integer*2 dist2(nres2,nres2)
c
	integer i,j
	real x,realdistance
c
c	for safety !
c
	do i=1,nres2
		do j=1,nres2
			dist2(i,j)=160
		end do
	end do
        do i=1,nres2
c
c		use bins 0.1 A 0-10; 0.4 A 10-20; 1 A 20-95
c		!minimal distance=1
c
                dist2(i,i)=1
                do j=i+1,nres2
                  x=realdistance(ca(1,i),ca(2,i),ca(3,i),
     $                          ca(1,j),ca(2,j),ca(3,j))
		  if(x.le.10.0) then
			dist2(i,j)=nint(x/0.1)
		  else if(x.le.20.0) then
			dist2(i,j)=100+nint((x-10.0)/0.4)
		  else if(x.lt.55.0) then
			dist2(i,j)=125+nint(x-20.0)
		  end if
                  dist2(j,i)=dist2(i,j)
                end do
        end do

        return
        end
c
c------------------------------------------------------------------------------
c
	subroutine getemax(e,emax,lclosed,iseg,jseg,xseg,yseg,seed)
	implicit none
	include 'parsizes.for'
	integer e,emax
	real ran
	logical lclosed(maxseg)
	integer iseg,jseg,xseg,yseg
	integer seed

	if(e.gt.emax) then
	    if(lclosed(iseg).and.(.not.lclosed(jseg))) then
		emax=e
		xseg=jseg
		yseg=iseg
	    else if(lclosed(jseg).and.(.not.lclosed(iseg))) then
		emax=e
		xseg=iseg
		yseg=jseg
	    else if(.not.lclosed(iseg).and.(.not.lclosed(jseg)))then
	        emax=e
		if(ran(seed).gt.0.5) then
			xseg=iseg
			yseg=jseg
		else
			xseg=jseg
			yseg=iseg
		end if
	    end if
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	function re_estimate(est,ns,seglist,ess,xseg,ex,ni,ci,start,mi,
     $ trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer re_estimate,trans(maxres0,maxseg)
	integer est,ns,seglist(maxseg),ess(maxseg,maxseg),xseg,mi(maxseg)
	integer ex(exdim),ni(maxseg),ci(maxres0,maxseg),
     $ start(maxseg,maxseg)
	logical lseqtl
c
	integer f,is,iseg,e
c
	f=est
	do is=1,ns
		f=f-ess(seglist(is),xseg)
	end do
	do is=1,ns
	  iseg=seglist(is)
	  call get_estimate(e,xseg,iseg,ni,ci,ex,start,mi,trans,lseqtl)
	  if(iseg.ne.xseg) then
		f=f+e+e
	  else
		f=f+e
	  end if
	end do

	re_estimate=f
c	write(*,*) 're-estimate:',f,est,xseg

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine bitpack(s,l,size)
	implicit none
	integer s,size
	logical l(31)
c
c	in: l
c	out: s, where l(i) -> ith bit in s
c
	integer a,i
c
	a=1
	s=0
	do i=1,size
		if(l(i)) s=s+a
		a=ishft(a,1)
	end do

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine bitunpack(s,l,b0,size)
	implicit none
	integer s,b0,size
	logical l(31)
c
c	in: s
c	out: l, where l(i) <- ith bit in s
c
	integer i,b
c
	b=b0
	do i=size,1,-1
		l(i)=(s.gt.b)
		s=iand(s,b)
		b=ishft(b,-1)
	end do

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine parsireadproteindata(nres,nseg,na,nb,secstr,
     $  segmentrange,
     $		checkrange,checkx,ca,ndom,node_type,node_child,domns,
     $		domseglist,iunit)
	implicit none
	include 'parsizes.for'
        integer ndom,node_child(2,maxdom),nres,nseg,iunit
        character node_type(maxdom),secstr(maxseg)
	integer domns(maxdom),domseglist(maxseg,maxdom),na,nb
	real ca(3,maxres)
	integer segmentrange(2,maxseg),checkrange(2,maxseg),checkx(maxseg)
c
	integer i,j,idom,iseg
c
	read(iunit,500) nres,nseg,na,nb,(secstr(i),i=1,nseg)
	do iseg=1,nseg
	  read(iunit,510) i,(segmentrange(j,i),j=1,2),
     $ (checkrange(j,i),j=1,2),checkx(i)
	end do
	read(iunit,520) ((ca(j,i),j=1,3),i=1,nres)
	read(iunit,500) ndom
	do idom=1,ndom
	  read(iunit,530) i,node_type(i),(node_child(j,i),j=1,2),domns(i),
     $		(domseglist(j,i),j=1,domns(i))
	end do

c	500	format(10x,4i5,2x,<maxseg>a1)

500	format(10x,4i5,2x,200a1)
510	format(6i10)
520	format(10f8.1)
c	530 format(i4,1x,a1,1x,3i4,<maxseg>i4)

530	format(i4,1x,a1,1x,3i4,200i4)

	return
	end
c
c-------------------------------------------------------------------------------
c
        function parsidistance(a1,a2,a3,b1,b2,b3)
        implicit none
        real a1,a2,a3,b1,b2,b3
        integer*2 parsidistance

        parsidistance=nint(10.0*sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+
     $          (a3-b3)*(a3-b3)))

        return
        end
c
c-----------------------------------------------------------------------
c
        function realdistance(a1,a2,a3,b1,b2,b3)
        implicit none
        real a1,a2,a3,b1,b2,b3
        real realdistance

        realdistance=sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+
     $          (a3-b3)*(a3-b3))

        return
        end
c
c-----------------------------------------------------------------------
c
        function distanceint4(a1,a2,a3,b1,b2,b3)
        implicit none
        real a1,a2,a3,b1,b2,b3
        integer distanceint4

	distanceint4=max(100,nint(1000.0*sqrt((a1-b1)*(a1-b1)
     $		+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))))
        return
        end
c
c-----------------------------------------------------------------------
c
	subroutine getlfix(idom,lfix,fixcutoff,ns,seglist,sf)
	implicit none
	include 'parsizes.for'
	integer idom
	logical lfix(maxseg,maxdom)
	real sf(maxseg,maxseg),fixcutoff
	integer ns,seglist(maxseg)
c
	integer is,iseg,js
	real x(maxseg)
c
	if(idom.eq.0) return
	do is=1,ns
		iseg=seglist(is)
		x(iseg)=0.0
		do js=1,ns
			if(js.ne.is) x(iseg)=x(iseg)+sf(iseg,seglist(js))
		end do
		lfix(iseg,idom)=(x(iseg).ge.fixcutoff)
	end do
c	if(ns.gt.1) write(*,500) idom,(seglist(is),x(seglist(is)),is=1,ns)
c
500	format('dom: ',i3,9(i3,f5.2),8(/10(i3,f5.2)))
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine segmentpairs(ali1,seglist,ns,idom,segment2,nseg2,
     $		checkrange,string,e,secstr,secstr2,nres2,checkx)
	implicit none
	include 'parsizes.for'
	integer ali1(maxseg),seglist(maxseg),ns,idom,segment2(maxres2)
	integer nseg2,checkrange(2,maxseg),e,nres2,checkx(maxseg)
	character*10 string
	character secstr(maxseg),secstr2(maxseg)
c
	integer npair,pair(2,90),i,k,x(0:maxseg),ires,jres,iseg,jseg,is
	logical lm
c
	npair=0
	do is=1,ns
		iseg=seglist(is)
		if(ali1(iseg).eq.nul) goto 19
		do jseg=0,nseg2
			x(jseg)=0
		end do
		do ires=checkrange(1,iseg),checkrange(2,iseg)
			jres=ali1(iseg)-checkx(iseg)+ires-checkrange(1,iseg)
			if((jres.ge.1).and.(jres.le.nres2)) then
				jseg=segment2(jres)
				x(jseg)=x(jseg)+1
			end if
		end do
		k=2
		lm=.false.
		if(secstr(iseg).eq.'H') k=6
		do jseg=1,nseg2
		  if (secstr2(jseg).eq.secstr(iseg)) then
		    if(x(jseg).ge.k) then
			npair=npair+1
			pair(1,npair)=iseg
			pair(2,npair)=jseg
			lm=.true.
		    end if
		  end if
		end do
		if(.not.lm) then
			npair=npair+1
			pair(1,npair)=iseg
			pair(2,npair)=0
		end if
19	end do
	!write(89,500) string,e,idom,npair,(pair(1,i),pair(2,i),i=1,npair)
c
500	format(a10,i10,2i5,5(2i4),10(/8(2i4)))
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine treehack(ndom,node_child,domns,domseglist)
	implicit none
	include 'parsizes.for'
	integer ndom
	integer node_child(2,maxdom),domns(maxdom),
     $ domseglist(maxseg,maxdom)
c
	integer i,j,idom,ndom0
c
c	separate multi-segment bottom nodes
c
c
c	write(*,*) 'original tree:'
c	do idom=1,ndom
c		write(*,500) idom,(node_child(j,idom),j=1,2),domns(idom),
c      $			(domseglist(j,idom),j=1,domns(idom))
c	end do
	idom=0
	do while(idom.lt.ndom)
	  idom=idom+1
	  if(domns(idom).gt.1.and.node_child(1,idom).eq.0) then
c         write(*,*) 'treehack:',ndom,idom,(domseglist(i,idom),
c     $ i=1,domns(idom))
		ndom=ndom+1
		node_child(1,idom)=ndom
		domns(ndom)=1
		domseglist(1,ndom)=domseglist(1,idom)
		do j=1,2
			node_child(j,ndom)=0
		end do
		ndom=ndom+1
		node_child(2,idom)=ndom
		domns(ndom)=domns(idom)-1
		do j=2,domns(idom)
			domseglist(j-1,ndom)=domseglist(j,idom)
		end do
		do j=1,2
			node_child(j,ndom)=0
		end do
	  end if
	end do
c	write(*,*) 'checked tree:'
c	do idom=1,ndom
c		write(*,500) idom,(node_child(j,idom),j=1,2),domns(idom),
c     $			(domseglist(j,idom),j=1,domns(idom))
c	end do
500	format(100i4)
	return
	end
c
c----------------------------------------------------------------------
c
