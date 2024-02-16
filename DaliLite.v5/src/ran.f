c This module/program is part of DaliLite (c) L. Holm 1999
c

	function ran(seed) 
	implicit none
	real ran 
	integer*4 seed,x,y,a,c
	data a,c,y/155,138,65536/

	x=mod (a*seed+c, y)
	seed=x
	ran=float(x)/float(y)

	return
	end
