DSET ^bvortex.dat
TITLE initialization t=1: bgd, t=2: truth, t>2: bgd+rnd 
UNDEF    -99999
XDEF  114 linear 0  2.5
YDEF  17  linear 0  2.5
ZDEF  1   linear 0  1
TDEF  33  linear 00z01apr85 1hr
vars 3 
u 0 99 truth
v 0 99 obs
z 0 99 height
endvars

