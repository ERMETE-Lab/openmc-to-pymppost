Cs-137 pencil beam on bare 3x3 inch NaI(Tl) - code-to-code verification
c ==================== cell cards ====================
c NaI(Tl) crystal: r = 3.81 cm, z = 0 to 7.62 cm, no housing / reflector / PMT
11   99  -3.67    -1  2  -3          imp:p=1     $ NaI(Tl) crystal (detector)
21    0           (1: -2: 3)  -4     imp:p=1     $ vacuum around the crystal
3     0            4                 imp:p=0     $ outside world, killed

c ==================== surface cards ====================
1    cz  3.81                        $ crystal radius (7.62 cm diameter)
2    pz  0.0                         $ crystal front face
3    pz  7.62                        $ crystal back face
4    so  50                          $ outer terminating sphere, r = 50 cm

c ==================== data cards ====================
mode  p
c --- NaI(Tl), 0.1 mol% Tl, weight fractions (negative), 3.67 g/cm3
m99  plib=04p
     11000  -0.153035
     53000  -0.845603
     81000  -0.001362
c --- Cs-137 (Ba-137m) gamma line, DDEP 2023: 661.657 keV
c     pencil beam on the crystal axis, starting 1 cm before the front face
sdef  pos=0 0 -1  cel=21  par=2  erg=0.661657  vec=0 0 1  dir=1
c --- photon physics: MCNPX defaults (coherent, Doppler broadening,
c     fluorescence, thick-target bremsstrahlung for electrons)
cut:p  j  0.001                      $ 1 keV photon cutoff, same as OpenMC
c --- MCNPX-PoliMi
c IPOL: (1)=0 source from SDEF; (2)-(4)=0 fission options, unused without
c       neutrons; (5)-(6) default; (7)=1 detector cell, (8)=its number
ipol  0 0 0 0 2j 1 11
c RPOL: (1) neutron, (2) photon energy-release thresholds [MeV]; 0 = write
c       every energy release, as the OpenMC collision file does
rpol  0 0
files 21 dumn1
prdmp 2j 1
nps   1e3
dbcn