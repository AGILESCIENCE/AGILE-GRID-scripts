#!/bin/bash
POSTFIX="-RUNTEST4"
#with GIFIXED
#jobfile="run_extended3.ll"
#without GI
jobfile="run_extended4.ll"

export SOURCE=HB21$POSTFIX L=88.752871 B=4.670431 R=1.19 A=''
sbatch $jobfile
sleep 1
export SOURCE=HB21$POSTFIX L=88.752871 B=4.670431 R=0.595 A=''
sbatch $jobfile
sleep 1
exit

export SOURCE=HESSJ1632-478$POSTFIX L=336.517018 B=0.121078 R=0.35 A=''
sbatch $jobfile
sleep 1

exit

export SOURCE=HESSJ1841-055$POSTFIX L=26.795483 B=-0.198336 R=0.62 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1841-055$POSTFIX L=25.867876 B=-0.362982 R=0.40 A=''
sbatch $jobfile
sleep 1

exit


export SOURCE=W28-HESSJ1800-240ABC$POSTFIX L=5.96 B=-0.39 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=W28-HESSJ1800-240ABC$POSTFIX L=5.96 B=-0.39 R=0.40 A=''
sbatch $jobfile
sleep 1
exit

export SOURCE=HESSJ1731-347$POSTFIX L=353.54 B=-0.68 R=0.34 A=''
sbatch $jobfile
sleep 1
exit

export SOURCE=RXJ1713.7-3946$POSTFIX L=347.335519 B=-0.472686 R=0.56 A=''
sbatch $jobfile
sleep 1
export SOURCE=RXJ1713.7-3946$POSTFIX L=347.34 B=-0.47 R=0.65 A=''
sbatch $jobfile
sleep 1
exit
export SOURCE=S147$POSTFIX L=180.244256 B=-1.499817 R=1.5 A=''
sbatch $jobfile
sleep 1
export SOURCE=S147$POSTFIX L=180.244256 B=-1.499817 R=0.75 A=''
sbatch $jobfile
sleep 1
exit


export SOURCE=W30$POSTFIX L=8.603773 B=-0.210545 R=0.37 A=''
sbatch $jobfile
sleep 1
export SOURCE=W30$POSTFIX L=8.40 B=-0.03 R=0.274 A=''
sbatch $jobfile
sleep 1
exit


export SOURCE=W44$POSTFIX L=34.653359 B=-0.388004 R=0.30 A=''
sbatch $jobfile
sleep 1

export SOURCE=W44$POSTFIX L=34.653359 B=-0.388004 R=0.15 A='' 
sbatch  $jobfile
sleep 1

export SOURCE=IC443$POSTFIX L=189.047650 B=3.033451 R=0.27 A=''
#sbatch $jobfile
sleep 1
export SOURCE=IC443$POSTFIX L=189.072910 B=2.917722 R=0.16 A=''
#sbatch $jobfile
sleep 1

exit

export SOURCE=VELAPS$POSTFIX L=263.58 B=-2.84 R=0.1 A=''
sbatch --partition=large $jobfile
sleep 1

exit

export SOURCE=HESSJ1825-137$POSTFIX L=17.566951 B=-0.453199 R=0.70 A=''
sbatch $jobfile
sleep 1

exit
exit
export SOURCE=CenALobes$POSTFIX L=309.167546 B=18.977150 R=2.5 A=''
sbatch $jobfile
sleep 1
export SOURCE=CenALobes$POSTFIX L=309.167546 B=18.977150  R=1.25 A=''
sbatch $jobfile
sleep 1


export SOURCE=CygnusCocoon$POSTFIX L=79.600721 B=1.396268 R=3.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=CygnusCocoon$POSTFIX L=80.953362 B=1.797498 R=1.8 A=''
sbatch $jobfile
sleep 1
export SOURCE=CygnusLoop$POSTFIX L=73.984182 B=-8.562410 R=3.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=CygnusLoop$POSTFIX L=73.980438 B=-10.570500 R=10 A=''
sbatch $jobfile
sleep 1
export SOURCE=HB21$POSTFIX L=88.752871 B=4.670431 R=1.19 A=''
sbatch $jobfile
sleep 1
export SOURCE=HB21$POSTFIX L=88.752871 B=4.670431 R=0.595 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1303-631$POSTFIX L=304.234773 B=-0.357503 R=0.24 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1303-631$POSTFIX L=304.213244 B=-0.334022 R=0.194 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1614-518$POSTFIX L=331.659006 B=-0.659082 R=0.42 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1614-518$POSTFIX L=331.519515 B=-0.581325 R=0.23 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1616-508$POSTFIX L=332.365145 B=-0.130941 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1616-508$POSTFIX L=332.390248 B=-0.141240 R=0.183 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1632-478$POSTFIX L=336.517018 B=0.121078 R=0.35 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1632-478$POSTFIX L=336.384259 B=0.190197 R=0.21 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1825-137$POSTFIX L=17.566951 B=-0.453199 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1825-137$POSTFIX L=17.710654 B=-0.696871 R=0.13 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1837-069$POSTFIX L=25.173028 B=-0.106905 R=0.33 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1837-069$POSTFIX L=25.177585 B=-0.115723 R=0.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1841-055$POSTFIX L=26.795483 B=-0.198336 R=0.62 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1841-055$POSTFIX L=25.867876 B=-0.362982 R=0.40 A=''
sbatch $jobfile
sleep 1
export SOURCE=LMC$POSTFIX L=278.842984 B=-32.850263 R=1.87 A=''
sbatch $jobfile
sleep 1
export SOURCE=LMC$POSTFIX L=279.55 B=-31.75 R=0.14 A=''
sbatch $jobfile
sleep 1
export SOURCE=MSH15-52$POSTFIX L=320.269499 B=-1.271354 R=0.25 A=''
sbatch $jobfile
sleep 1
export SOURCE=MSH15-52$POSTFIX L=320.33 B=-1.19 R=0.11 A=''
sbatch $jobfile
sleep 1
export SOURCE=PuppisA$POSTFIX L=260.317182 B=-3.276471 R=0.37 A=''
sbatch $jobfile
sleep 1
export SOURCE=PuppisA$POSTFIX L=260.317182 B=-3.276471 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=RXJ1713.7-3946$POSTFIX L=347.335519 B=-0.472686 R=0.56 A=''
sbatch $jobfile
sleep 1
export SOURCE=RXJ1713.7-3946$POSTFIX L=347.34 B=-0.47 R=0.65 A=''
sbatch $jobfile
sleep 1
export SOURCE=S147$POSTFIX L=180.244256 B=-1.499817 R=1.5 A=''
sbatch $jobfile
sleep 1
export SOURCE=S147$POSTFIX L=180.244256 B=-1.499817 R=0.75 A=''
sbatch $jobfile
sleep 1
export SOURCE=SMC$POSTFIX L=302.144948 B=-44.416694 R=1.35 A=''
sbatch $jobfile
sleep 1
export SOURCE=SMC$POSTFIX L=302.144948 B=-44.416694 R=0.67 A=''
sbatch $jobfile
sleep 1
export SOURCE=VELAJr$POSTFIX L=266.490859  B=-1.233154 R=1.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=VELAJr$POSTFIX L=266.28  B=-1.24 R=1.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=VELAX$POSTFIX L=263.33208 B=-3.10605 R=0.91 A=''
sbatch $jobfile
sleep 1
export SOURCE=VELAX$POSTFIX L=263.86 B=-3.09 R=0.48 A=''
sbatch $jobfile
sleep 1
export SOURCE=W51C$POSTFIX L=49.115799 B=-0.461602 R=0.375 A=''
sbatch $jobfile
sleep 1
export SOURCE=W51C$POSTFIX L=49.12 B=-0.36 R=0.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=gammaCygni$POSTFIX L=78.240815 B=2.196722 R=0.63 A=''
sbatch $jobfile
sleep 1
export SOURCE=gammaCygni$POSTFIX L=78.33 B=2.49 R=0.23 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTA1$POSTFIX L=119.6 B=10.4 R=0.30 A=''
sbatch $jobfile
sleep 1
export SOURCE=Geminga$POSTFIX L=195.33 B=3.77 R=1.3 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1018-589B$POSTFIX L=284.11 B=-1.9 R=0.15 A=''
sbatch $jobfile
sleep 1
export SOURCE=Westerlund2$POSTFIX L=284.21 B=-0.41 R=0.18 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1026-582$POSTFIX L=284.79 B=-0.53 R=0.14 A=''
sbatch $jobfile
sleep 1
export SOURCE=SNRG292.2-00.5$POSTFIX L=292.1 B=-0.49 R=0.1 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1356-645$POSTFIX L=309.81 B=-2.5 R=0.2 A=''
sbatch $jobfile
sleep 1
export SOURCE=Kookaburra-Rabbit$POSTFIX L=313.24 B=0.14 R=0.082 A=''
sbatch $jobfile
sleep 1
export SOURCE=Kookaburra-Pulsar$POSTFIX L=313.55 B=0.26 R=0.055 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1427-608$POSTFIX L=314.4 B=-0.15 R=0.04 A=''
sbatch $jobfile
sleep 1
export SOURCE=RCW86$POSTFIX L=315.41 B=-2.31 R=0.41 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1457-593$POSTFIX L=318.36 B=-0.44 R=0.31 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1458-608$POSTFIX L=317.74 B=-1.71 R=0.17 A=''
sbatch $jobfile
sleep 1
export SOURCE=SN1006-SW$POSTFIX L=327.86 B=15.34 R=0.13 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1503-582$POSTFIX L=319.61 B=0.29 R=0.26 A=''
sbatch $jobfile
sleep 1
export SOURCE=SN1006-NE$POSTFIX L=327.84 B=14.56 R=0.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1507-622$POSTFIX L=317.94 B=-3.5 R=0.15 A=''
sbatch $jobfile
sleep 1
export SOURCE=G327.1-1.1$POSTFIX L=327.15 B=-1.08 R=0.03 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1626-490$POSTFIX L=334.77 B=0.04 R=0.1 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1634-472$POSTFIX L=337.1 B=0.21 R=0.11 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1640-465$POSTFIX L=338.31 B=-0.03 R=0.045 A=''
sbatch $jobfile
sleep 1
export SOURCE=Westerlund1$POSTFIX L=339.54 B=-0.36 R=1.1 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1702-420$POSTFIX L=344.3 B=-0.19 R=0.3 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1708-410$POSTFIX L=345.68 B=-0.48 R=0.08 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1708-443$POSTFIX L=343.05 B=-2.38 R=0.29 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTB37B$POSTFIX L=348.63 B=0.38 R=0.06 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTB37A$POSTFIX L=348.38 B=0.1 R=0.072 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1718-385$POSTFIX L=348.83 B=-0.49 R=0.15 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1729-345$POSTFIX L=353.44 B=-0.13 R=0.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1731-347$POSTFIX L=353.54 B=-0.68 R=0.27 A=''
sbatch $jobfile
sleep 1
export SOURCE=GalicticCenterRidge$POSTFIX L=359.94 B=-0.05 R=2.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1745-303$POSTFIX L=358.7 B=-0.65 R=0.21 A=''
sbatch $jobfile
sleep 1
export SOURCE=Terzan5$POSTFIX L=3.78 B=1.72 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=W28-HESSJ1800-240ABC$POSTFIX L=5.96 B=-0.39 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1804-216$POSTFIX L=8.35 B=-0.01 R=0.274 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1808-204$POSTFIX L=9.95 B=-0.25 R=0.138  A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1809-193$POSTFIX L=11.18 B=-0.09 R=0.53 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1813-178$POSTFIX L=12.81 B=-0.03 R=0.036 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1831-098$POSTFIX L=21.85 B=-0.11 R=0.15 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1834-087$POSTFIX L=23.24 B=-0.33 R=0.09 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1843-033$POSTFIX L=29.03 B=0.36 R=1.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1848-018$POSTFIX L=31.0 B=-0.17 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=IGRJ18490-0000$POSTFIX L=32.63 B=0.52 R=1.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1857+026$POSTFIX L=36 B=-0.07 R=0.20 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1858+020$POSTFIX L=35.57 B=-0.59 R=0.08 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ1908+06$POSTFIX L=40.28 B=-0.69 R=0.44 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1912+101$POSTFIX L=44.39 B=-0.08 R=0.27 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ2019+37$POSTFIX L=74.82 B=0.41 R=0.75 A=''
sbatch $jobfile
sleep 1
export SOURCE=VERJ2019+368$POSTFIX L=75.04 B=0.28 R=0.34 A=''
sbatch $jobfile
sleep 1
export SOURCE=VERJ2019+407$POSTFIX L=78.33 B=2.49 R=0.23 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ2031+41$POSTFIX L=79.53 B=0.63 R=1.8 A=''
sbatch $jobfile
sleep 1
export SOURCE=TeVJ2032+4130$POSTFIX L=80.24 B=1.17 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=G106.3+2.7$POSTFIX L=106.34 B=2.71 R=0.27 A=''
sbatch $jobfile
sleep 1
export SOURCE=Boomerang$POSTFIX L=106.57 B=2.91 R=1.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTA1$POSTFIX L=119.6 B=10.4 R=0.60 A=''
sbatch $jobfile
sleep 1
export SOURCE=Geminga$POSTFIX L=195.33 B=3.77 R=2.6 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1018-589B$POSTFIX L=284.11 B=-1.9 R=0.30 A=''
sbatch $jobfile
sleep 1
export SOURCE=Westerlund2$POSTFIX L=284.21 B=-0.41 R=0.36 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1026-582$POSTFIX L=284.79 B=-0.53 R=0.28 A=''
sbatch $jobfile
sleep 1
export SOURCE=SNRG292.2-00.5$POSTFIX L=292.1 B=-0.49 R=0.2 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1356-645$POSTFIX L=309.81 B=-2.5 R=0.4 A=''
sbatch $jobfile
sleep 1
export SOURCE=Kookaburra-Rabbit$POSTFIX L=313.24 B=0.14 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=Kookaburra-Pulsar$POSTFIX L=313.55 B=0.26 R=0.11 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1427-608$POSTFIX L=314.4 B=-0.15 R=0.08 A=''
sbatch $jobfile
sleep 1
export SOURCE=RCW86$POSTFIX L=315.41 B=-2.31 R=0.82 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1457-593$POSTFIX L=318.36 B=-0.44 R=0.62 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1458-608$POSTFIX L=317.74 B=-1.71 R=0.34 A=''
sbatch $jobfile
sleep 1
export SOURCE=SN1006-SW$POSTFIX L=327.86 B=15.34 R=0.26 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1503-582$POSTFIX L=319.61 B=0.29 R=0.52 A=''
sbatch $jobfile
sleep 1
export SOURCE=SN1006-NE$POSTFIX L=327.84 B=14.56 R=0.24 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1507-622$POSTFIX L=317.94 B=-3.5 R=0.30 A=''
sbatch $jobfile
sleep 1
export SOURCE=G327.1-1.1$POSTFIX L=327.15 B=-1.08 R=0.06 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1626-490$POSTFIX L=334.77 B=0.04 R=0.2 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1634-472$POSTFIX L=337.1 B=0.21 R=0.22 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1640-465$POSTFIX L=338.31 B=-0.03 R=0.01 A=''
sbatch $jobfile
sleep 1
export SOURCE=Westerlund1$POSTFIX L=339.54 B=-0.36 R=2.2 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1702-420$POSTFIX L=344.3 B=-0.19 R=0.6 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1708-410$POSTFIX L=345.68 B=-0.48 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1708-443$POSTFIX L=343.05 B=-2.38 R=0.60 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTB37B$POSTFIX L=348.63 B=0.38 R=0.12 A=''
sbatch $jobfile
sleep 1
export SOURCE=CTB37A$POSTFIX L=348.38 B=0.1 R=0.14 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1718-385$POSTFIX L=348.83 B=-0.49 R=0.30 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1729-345$POSTFIX L=353.44 B=-0.13 R=0.24 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1731-347$POSTFIX L=353.54 B=-0.68 R=0.34 A=''
sbatch $jobfile
sleep 1
export SOURCE=GalicticCenterRidge$POSTFIX L=359.94 B=-0.05 R=4.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1745-303$POSTFIX L=358.7 B=-0.65 R=0.42 A=''
sbatch $jobfile
sleep 1
export SOURCE=Terzan5$POSTFIX L=3.78 B=1.72 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=W28-HESSJ1800-240ABC$POSTFIX L=5.96 B=-0.39 R=0.64 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1804-216$POSTFIX L=8.35 B=-0.01 R=0.348 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1808-204$POSTFIX L=9.95 B=-0.25 R=0.28  A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1809-193$POSTFIX L=11.18 B=-0.09 R=1.06 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1813-178$POSTFIX L=12.81 B=-0.03 R=0.080 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1831-098$POSTFIX L=21.85 B=-0.11 R=0.30 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1834-087$POSTFIX L=23.24 B=-0.33 R=0.18 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1843-033$POSTFIX L=29.03 B=0.36 R=2.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1848-018$POSTFIX L=31.0 B=-0.17 R=0.64 A=''
sbatch $jobfile
sleep 1
export SOURCE=IGRJ18490-0000$POSTFIX L=32.63 B=0.52 R=2.0 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1857+026$POSTFIX L=36 B=-0.07 R=0.40 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1858+020$POSTFIX L=35.57 B=-0.59 R=0.16 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ1908+06$POSTFIX L=40.28 B=-0.69 R=0.88 A=''
sbatch $jobfile
sleep 1
export SOURCE=HESSJ1912+101$POSTFIX L=44.39 B=-0.08 R=0.52 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ2019+37$POSTFIX L=74.82 B=0.41 R=1.5 A=''
sbatch $jobfile
sleep 1
export SOURCE=VERJ2019+368$POSTFIX L=75.04 B=0.28 R=0.68 A=''
sbatch $jobfile
sleep 1
export SOURCE=VERJ2019+407$POSTFIX L=78.33 B=2.49 R=0.46 A=''
sbatch $jobfile
sleep 1
export SOURCE=MGROJ2031+41$POSTFIX L=79.53 B=0.63 R=3.6 A=''
sbatch $jobfile
sleep 1
export SOURCE=TeVJ2032+4130$POSTFIX L=80.24 B=1.17 R=0.32 A=''
sbatch $jobfile
sleep 1
export SOURCE=G106.3+2.7$POSTFIX L=106.34 B=2.71 R=0.52 A=''
sbatch $jobfile
sleep 1
export SOURCE=Boomerang$POSTFIX L=106.57 B=2.91 R=2.0 A=''
sbatch $jobfile
sleep 1
