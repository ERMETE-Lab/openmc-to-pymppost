stilbene detector - DD source
c ==================== cell cards ====================
11   99  -1.16    -1  2  -3          imp:n,p=1     $ stilbene
21    0            (1: -2: 3)  -5    imp:n,p=1     $ void, interno alla sfera fittizia
22    0            5  -4              imp:n,p=1     $ void, tra sfera fittizia e sfera esterna
3     0            4                 imp:n,p=0     $ outside world, killed

c ==================== surface cards ====================
1    cz  2.54                        $ detector radius (5.08 cm diameter)
2    pz  0.0                         $ detector bottom
3    pz  5.08                        $ detector top
5    so  15                          $ sfera fittizia per tally di corrente, r = 15 cm
4    so  50                          $ outer terminating sphere, r = 50 cm

c ==================== data cards ====================
c --- trans-stilbene
m99  nlib=71c plib=04p pnlib=70u
     1001  .4615
     6000  .5385
c --- nuclidi "puri" per i moltiplicatori FM (scattering per isotopo)
m101 nlib=71c
     1001  1
m102 nlib=71c
     6000  1
c
c --- cf-252 point source, 10 cm below detector, monodirectional
sdef  pos=0 0 -10  cel=21  par=1  erg=2.45  vec=0 0 1  dir=1
c
fc14 cell neutron flux in detector
f14:n  11
e14  0  0.01  199i  2.5
c
fc24 elastic scattering (MT=2) reaction rate on H-1 in detector
f24:n  11
fm24   0.046507  101  2
e24  0  0.01  199i  2.5
c
fc34 elastic scattering (MT=2) reaction rate on C in detector
f34:n  11
fm34   0.054258  102  2
e34  0  0.01  199i  2.5
c
fc41 neutron current across fictitious sphere (r=15 cm)
f41:n  5
e41  0  0.01  199i  2.5
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