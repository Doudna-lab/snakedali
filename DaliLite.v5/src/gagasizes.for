c
c	include file for gaga*
c
	integer popsize,maxres,maxpair,maxd,maxres0,maxdom
	parameter(popsize=1)
	parameter(maxd=2000)
	parameter(maxpair=6000)
	parameter(maxres=5000,maxres0=maxres,maxdom=210)

c
        character machine
        real enveloperadius,tripletscorecutoff,killcutoff(3),d0,pblock
        real fokcutoff,putcutoff,seedcutoff,refitol,hack0,hack1
        integer plimit,p600,p720,iounit1,geneblocksize,outblocksize
        integer corecutoff,itrim(3),ikill(3),nclone(3),maxstep,nrepeat
        integer printpopsize,nbest
        logical ltop,lpara,lali,lpretty,lsavegenes,lcalpha,lela
        logical lplottriplets
        integer*4 seed
c
c
c
	real weight(0:100)
	logical ldebug,lverb,lprintparams,ltube
	common /param/lverb,weight,d0,enveloperadius,
     $  tripletscorecutoff,killcutoff,pblock,
     $  corecutoff,fokcutoff,putcutoff,plimit,p600,p720,iounit1,
     $  geneblocksize,outblocksize,maxstep,itrim,ikill,nclone,
     $  ltop,lpara,lali,lpretty,lsavegenes,lcalpha,seed,nrepeat,
     $	nbest,printpopsize,seedcutoff,lela,refitol,hack0,hack1,
     $	lplottriplets,ldebug,lprintparams,ltube,
     $  machine

	integer dim1,dim23,dim4,dim5,dim6
	parameter(dim1=maxres*(maxres+maxres+1))
	parameter(dim23=dim1)
        parameter(dim4=maxres*maxres)
	parameter(dim5=maxres*maxres)
	parameter(dim6=maxres*maxres)
