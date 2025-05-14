Inc 69.0
End date string: 2025-03-25T00:00:00.000Z
Phase angle target: 180.0°
TEST
Loaded time/orbit data from:
 - maneuverOrderFile => DATE: 2025-04-30T21:04:18.672Z
 - commandDataFile => SMA: 1.5E7 m, ECC: 0.01, etc.
 - ergolsFile => DryMass: 2478.0 kg, Ergol: 3550.357775
 - maneuvFile => manoeuverRelativeDate: -3627057.2793278694
 - timeIntermediateParametersFile => ManeuverType: Inclinaison, PhaseAngle: 180.0°
Using ergols input file: Ergols.txt
Using ergols consumption output file: ConsommationErgols.txt
Mode parameter: blue1
isMassCalculation: false
CommandData.txt
{AoP=98.0, ECC=0.01, SMA=15000.0, ManeuverType=Phasage, DELTA_THETA=180, ISP=300, Ergol mass=300, DATE=, Dry Mass=800, MeanAnom=126.0, Thurst=0, RAAN=45.0, INC=69.0, SurfRef=0}
Using ergols input file: Ergols.txt
Masse à vide : 2478.0 kg
Masse d'ergols : 3550.357775 kg
Impulsion spécifique (ISP) : 3000.0 s
SMA initial : 15000.0 km

=== PHASING MANOEUVRE  ==================================================
Dry mass              : 2 478,000 kg
Propellant on‑board   : 3 550,358 kg
ISP                   : 3 000,0 s
Initial wet mass m0   : 6 028,358 kg
Reference apsis date  : 2025-03-19T19:04:18.67286182652058Z
Burn‑1 epoch (t0)     : 2025-03-19T21:33:21.39267213058472Z

[state@apsis]
  Epoch   : 2025-03-19T19:04:18.67286182652058Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : -262,000 °
  M       : 126,000 °
  ν (TA)  : 126,920 °

[state@initial]
  Epoch   : 2025-03-19T21:33:21.39267213058472Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 98,000 °
  M       : 302,086 °
  ν (TA)  : 301,108 °
Δθ commanded          : +180,0000 deg
Natural drift         : +176,0858 deg
Δθ still to perform   : -180,0000 deg
k_int                 : 1
T_ph (initial guess)  : 457,075 min
Phasing orbit (first) : a=19655,560 km  e=0,236857  rp=15000,000 km  ra=24311,121 km
ΔV1 / ΔV2             : +578,079 / -578,079 m/s
Propellant expected   : 232,310 kg

[state@afterBurn1]
  Epoch   : 2025-03-19T21:33:21.39467213058472Z
  Mass    : 5911,061 kg
  a       : 19685,954 km
  e       : 0,24216879
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 41,612 °
  M       : 358,518 °
  ν (TA)  : 357,496 °
Newton refine  it=0  Δθ=-1,407097 deg  dt=-71,5 s  → Tph=458,266 min
Newton refine  it=1  Δθ=+0,173064 deg  dt=+8,8 s  → Tph=458,120 min
Newton refine  it=2  Δθ=-0,021332 deg  dt=-1,1 s  → Tph=458,138 min
Newton refine  it=3  Δθ=+0,002629 deg  dt=+0,1 s  → Tph=458,136 min
Newton refine  it=4  Δθ=-0,000324 deg  dt=-0,0 s  → Tph=458,136 min
Newton refine  it=5  |Δθ|=0,000040 deg  ✓

[state@beforeBurn2]
  Epoch   : 2025-03-20T05:11:29.5566834083306Z
  Mass    : 5911,061 kg
  a       : 19685,954 km
  e       : 0,24216879
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 41,612 °
  M       : 358,518 °
  ν (TA)  : 357,496 °
SMA (Semi-major axis): 19685.953814849425 km

[Reference state (no maneuver)]
  Epoch   : 2025-03-20T05:11:29.5566834083306Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : -262,000 °
  M       : 123,339 °
  ν (TA)  : 124,290 °

[state@final]
  Epoch   : 2025-03-20T05:11:29.5586834083306Z
  Mass    : 5796,047 kg
  a       : 15000,000 km
  e       : 0,00999976
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 97,999 °
  M       : 302,087 °
  ν (TA)  : 301,109 °
Co-location check     : Δr = 21,341 m   Δθ = 0,000082 deg
Phase achieved vs ref : +178,7470 deg
SMA (Semi-major axis): 14999.999962131735 km