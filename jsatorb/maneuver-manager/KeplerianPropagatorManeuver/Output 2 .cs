/opt/jdk-23.0.1/bin/java -javaagent:/home/simulateurspatial/.local/share/JetBrains/Toolbox/apps/intellij-idea-ultimate/lib/idea_rt.jar=46041:/home/simulateurspatial/.local/share/JetBrains/Toolbox/apps/intellij-idea-ultimate/bin -Dfile.encoding=UTF-8 -Dsun.stdout.encoding=UTF-8 -Dsun.stderr.encoding=UTF-8 -classpath /home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/out/production/OrbitalManeuver:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/orekit-11.3.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-core-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-geometry-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-ode-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-fitting-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-optim-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-filtering-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/hipparchus-stat-2.3.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/json-20250107.jar:/home/simulateurspatial/jsatorb/maneuver-manager/KeplerianPropagatorManeuver/lib/org.eclipse.paho.client.mqttv3-1.2.5.jar OrbitalManeuver
Inc 69.0
End date string: 2025-03-25T00:00:00.000Z
180.0test
TEST
Loaded time/orbit data from:
 - maneuverOrderFile => DATE: 2025-04-30T21:04:18.672Z
 - commandDataFile => SMA: 1.5E7 m, ECC: 0.01, etc.
 - ergolsFile => DryMass: 2478.0 kg, Ergol: 3550.357775
 - maneuvFile => manoeuverRelativeDate: -3627057.2793278694
 - timeIntermediateParametersFile => ManeuverType: Inclinaison, PhaseAngle: 10313.240312354817°
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
Delta theta in degrees: 180.0 °

[State@apside]
  Epoch   : 2025-03-19T19:04:18.67286182652058Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : -262,000 °
  M       : 126,000 °
  ν (TA)  : 126,920 °

[State@initialDate]
  Epoch   : 2025-03-19T21:33:21.39267213058472Z
  Mass    : 1000,000 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 98,000 °
  M       : 302,086 °
  ν (TA)  : 301,108 °

[state0]
  Epoch   : 2025-03-19T21:33:21.39267213058472Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 98,000 °
  M       : 302,086 °
  ν (TA)  : 301,108 °
Delta theta in radians: 3.141592653589793 rad
Delta theta in degrees: 180.0 °
Kint: 2

[Reference state (no maneuver)]
  Epoch   : 2025-03-20T05:10:25.91855093184992Z
  Mass    : 6028,358 kg
  a       : 15000,000 km
  e       : 0,01000000
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : -262,000 °
  M       : 122,086 °
  ν (TA)  : 123,050 °
Warning: Swapped ra and rp to maintain physical constraint
Original orbit: r=15 000 km, period=304,72 min
Phasing orbit: a=19 656 km, rp=15 000 km, ra=24 311 km, period=457,08 min, e=0,23686
Delta-Vs: dv1=578,08 m/s, dv2=-578,08 m/s, total=1 156,16 m/s
Mass: initial=6 028,36 kg, after burn 1=5 911,06 kg, final=5 796,05 kg, fuel used=232,31 kg

[State@afterBurn1]
  Epoch   : 2025-03-19T21:33:21.39367213058472Z
  Mass    : 5911,061 kg
  a       : 19685,954 km
  e       : 0,24216879
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 41,612 °
  M       : 358,518 °
  ν (TA)  : 357,496 °

[State@beforeBurn2]
  Epoch   : 2025-03-20T05:10:25.91855093184992Z
  Mass    : 5911,061 kg
  a       : 19685,954 km
  e       : 0,24216879
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 41,612 °
  M       : 357,684 °
  ν (TA)  : 356,089 °
Payload for first burn: {enddate=2025-03-20T05:10:25.919Z, startdate=2025-03-19T21:33:21.393Z}
SMA (Semi-major axis): 19685.953814849418 km

[State@afterBurn1]
  Epoch   : 2025-03-19T21:33:21.39367213058472Z
  Mass    : 5911,061 kg
  a       : 19685,954 km
  e       : 0,24216879
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 41,612 °
  M       : 358,518 °
  ν (TA)  : 357,496 °

[State@final]
  Epoch   : 2025-03-20T05:10:25.91955093184993Z
  Mass    : 5796,047 kg
  a       : 15000,811 km
  e       : 0,01421924
  i       : 69,000 °
  RAAN    : 45,000 °
  ω (PA)  : 107,763 °
  M       : 291,465 °
  ν (TA)  : 289,939 °
Phase change: target=180,00°, achieved=169,38°, theoretical=-180,00°
Phase change: target=180,00°, achieved=169,38°, absolute=169,38°
Mean Anomaly: initial=302,09°, final=291,46°, reference=122,09°
PHASInG END 2025-03-20T05:10:25.91855093184992Z
Payload for second burn: {enddate=2025-03-25T00:00:00Z, startdate=2025-03-20T05:10:25.920Z}
SMA (Semi-major axis): 15000.810958381042 km
postManeuverDate 2025-03-20T05:10:25.91955093184993Z
Phasing maneuver complete. Execution time: 0,228 s

Process finished with exit code 0
