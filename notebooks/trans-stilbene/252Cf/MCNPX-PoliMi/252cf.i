stilbene detector
c ==================== cell cards ====================
11   99  -1.16    -1  2  -3          imp:n,p=1     $ stilbene
2     0            (1: -2: 3)  -4    imp:n,p=1     $ void
3     0            4                 imp:n,p=0     $ outside world, killed

c ==================== surface cards ====================
1    cz  2.54                        $ detector radius (5.08 cm diameter)
2    pz  0.0                         $ detector bottom
3    pz  5.08                        $ detector top
4    so  50                          $ outer terminating sphere, r = 50 cm

c ==================== data cards ====================
c --- trans-stilbene
m99  nlib=70c plib=04p pnlib=70u
     1001  .4615
     6000  .5385
c
c --- cf-252 point source, 10 cm below detector, monodirectional
sdef  pos=0 0 -10  cel=2  par=1  erg=d1  vec=0 0 1  dir=1
c watt fission spectrum approximation for cf-252
c (a=1.025 mev, b=2.926 mev^-1 - verify against reference)
sp1  -3  1.025  2.926
c
fc14 cell neutron flux in detector
f14:n  11
e14  0  0.01  199i  20
c
phys:n  j  20
phys:p  3j  -1
mode    n p
ipol  0 0 0 0 2j 1 11    $ one detector cell 11
rpol  0.01 100           $ threshold for not writing photon events
files 21 dumn1
prdmp 2j 1
nps   4e3
dbcn