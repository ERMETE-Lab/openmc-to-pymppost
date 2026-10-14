Eu152 pencil beam on bare 3x3 inch NaI(Tl) - code-to-code verification
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
c --- Eu-152 photon source: LNHB/DDEP recommended data (Vanin, de Castro,
c     Browne, 2004). Gamma lines >= 0.1 per 100 decays + Sm/Gd K X-rays.
c     48 lines, 2.3173 photons per decay.
c     Pencil beam along +z from 1 cm before the crystal front face.
sdef  pos=0 0 -1  cel=21  par=2  erg=d1  vec=0 0 1  dir=1
si1 L 0.1217817 0.2446974 0.2959387 0.3294250 0.3442785 0.3677891 0.4111165
     0.4160480 0.4439650 0.4886792 0.5034740 0.5639900 0.5664420 0.5862650
     0.6564890 0.6746750 0.6786230 0.6886700 0.7193490 0.7649000 0.7789045
     0.8104510 0.8415740 0.8673800 0.9193370 0.9263170 0.9633900 0.9640790
     1.0052720 1.0840000 1.0858370 1.0897370 1.1091740 1.1120760 1.2129480
     1.2499380 1.2927780 1.2991420 1.4080130 1.4576430 1.5281030 0.0395229
     0.0401186 0.0454130 0.0465750 0.0423093 0.0429967 0.0486970
sp1   2.841000e-01 7.550000e-02 4.420000e-03 1.290000e-03 2.659000e-01
     8.620000e-03 2.238000e-02 1.090000e-03 3.120000e-02 4.139000e-03
     1.533000e-03 4.570000e-03 1.310000e-03 4.620000e-03 1.437000e-03
     1.871000e-03 4.700000e-03 8.410000e-03 3.270000e-03 1.900000e-03
     1.297000e-01 3.170000e-03 1.630000e-03 4.243000e-02 4.290000e-03
     2.730000e-03 1.341000e-03 1.450000e-01 6.650000e-03 2.440000e-03
     1.013000e-01 1.730000e-02 1.860000e-03 1.341000e-01 1.416000e-02
     1.860000e-03 1.040000e-03 1.633000e-02 2.085000e-01 4.980000e-03
     2.810000e-03 2.080000e-01 3.770000e-01 1.178000e-01 3.040000e-02
     2.430000e-03 4.370000e-03 1.380000e-03
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