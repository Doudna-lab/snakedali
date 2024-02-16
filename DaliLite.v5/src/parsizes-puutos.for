	integer maxstack2,maxdom,maxres0,maxres1,maxres2,exdim
	
	parameter (maxdom=2060)
	parameter (maxres0=3000,maxres1=3000,maxres2=3000,exdim=350000)
	integer maxstack,maxseg,maxatm,maxres,maxgrid,maxrot
	parameter (maxstack=10000,maxstack2=10000)
	parameter (maxseg=320)
	integer maxprot,maxpair, maxnoccs
	integer weight(1000)
	integer scoretable(160,400)
	integer infinit,nul
	parameter (infinit=1e7,nul=-99)
	parameter (maxprot=1000000)
	parameter (maxatm=5000)
	parameter (maxres=12000)
	parameter (maxgrid=20)
	parameter (maxrot=500)
	parameter (maxpair=50000)
	parameter (maxnoccs=300000)
	integer rankcount,callcount,testcount

	integer boxdim
	parameter(boxdim=501)

	common weight,scoretable,rankcount,callcount,testcount
  
