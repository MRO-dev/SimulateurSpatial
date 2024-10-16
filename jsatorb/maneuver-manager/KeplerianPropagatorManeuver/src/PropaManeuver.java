//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.util.Locale;
import org.hipparchus.util.FastMath;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.MathUtils;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.LofOffset;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.forces.ForceModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.forces.maneuvers.ImpulseManeuver;
import org.orekit.forces.maneuvers.Maneuver;
import org.orekit.forces.maneuvers.propulsion.BasicConstantThrustPropulsionModel;
import org.orekit.forces.maneuvers.propulsion.ThrustPropulsionModel;
import org.orekit.forces.maneuvers.trigger.DateBasedManeuverTriggers;
import org.orekit.forces.maneuvers.trigger.ManeuverTriggers;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.LOFType;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

public class PropaManeuver {
    public PropaManeuver() {
    }

    public static void main(String[] args) throws NullPointerException, ClassCastException, IOException, ParseException {
        double start = (double)System.currentTimeMillis();

        try {
            int num = 1;
            String DATE = (String)Files.readAllLines(Paths.get("Data.txt")).get(num);
            num = 2;
            double SMA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 3;
            double ECC = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 4;
            double INC = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 5;
            double RAAN = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 6;
            double PA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 7;
            double ANO = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            double sma = SMA * 1000.0;
            double ecc = ECC;
            double inc = FastMath.toRadians(INC);
            double raan = FastMath.toRadians(RAAN);
            double aop = FastMath.toRadians(PA);
            double ano = FastMath.toRadians(ANO);
            double MU = 3.986004418E14;
            Frame GCRF = FramesFactory.getEME2000();
            num = 8;
            double DRYMASS = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            double dryMass = DRYMASS;
            num = 10;
            double THRUST = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 11;
            double ISP = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 12;
            double ERGOL = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 13;
            double DURA = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            num = 14;
            String maneuverType = (String)Files.readAllLines(Paths.get("Data.txt")).get(num);
            num = 15;
            String integratorType = (String)Files.readAllLines(Paths.get("Data.txt")).get(num);
            num = 16;
            double DV = Double.parseDouble((String)Files.readAllLines(Paths.get("Data.txt")).get(num));
            System.out.println("maneuverType" + maneuverType);
            System.out.println("integratorType" + integratorType);

            try {
                int item = 1;
                double manoeuverRelativeDate = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 2;
                double DVx = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 3;
                double DVy = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 4;
                double DVz = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                item = 5;
                double durationOfManoeuver = Double.parseDouble((String)Files.readAllLines(Paths.get("Maneuv.txt")).get(item));
                Vector3D direction = new Vector3D(DVx, DVy, DVz);
                System.out.println("Vecteur direction X= " + direction.getX() + " Y= " + direction.getY() + " Z= " + direction.getZ());

                try {
                    File home = new File("/app/maneuver-manager/");
                    File orekitData = new File(Paths.get("orekit-data").toString());
                    if (!orekitData.exists()) {
                        System.err.format(Locale.US, "Failed to find %s folder%n", orekitData.getAbsolutePath());
                        System.err.format(Locale.US, "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n", "orekit-data-master.zip", "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip", home.getAbsolutePath());
                        System.exit(1);
                    }

                    DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
                    manager.addProvider(new DirectoryCrawler(orekitData));
                    Frame eme2000 = FramesFactory.getEME2000();
                    AbsoluteDate dateTLE = new AbsoluteDate(DATE, TimeScalesFactory.getUTC());
                    System.out.println("DATE :" + DATE);
                    System.out.println("dateTLE :" + dateTLE);
                    Orbit iniOrbit = new KeplerianOrbit(sma, ecc, inc, aop, raan, ano, PositionAngle.MEAN, eme2000, dateTLE, 3.986004418E14);
                    double initMass = ERGOL + dryMass;
                    SpacecraftState initialState = new SpacecraftState(iniOrbit, initMass);
                    System.out.println("iniOrbit:" + iniOrbit);
                    System.out.println("initMass :" + initMass);
                    AbsoluteDate manoeuverStartDate = dateTLE.shiftedBy(manoeuverRelativeDate);
                    AbsoluteDate manoeuverEndDate = manoeuverStartDate.shiftedBy(durationOfManoeuver);
                    OrbitType orbitType = OrbitType.CIRCULAR;
                    if (ecc > 0.01) {
                        orbitType = OrbitType.KEPLERIAN;
                    } else {
                        orbitType = OrbitType.CIRCULAR;
                    }

                    NumericalPropagator propagator;
                    double stepSize;
                    if (integratorType.equals("DP")) {
                        System.out.println("Integrator : DormandPrince");
                        double[][] tol = NumericalPropagator.tolerances(1.0, iniOrbit, orbitType);
                        stepSize = 1.0;
                        double maxStep = stepSize * 100000.0;
                        double minStep = stepSize / 1000.0;
                        double initialStepSize = minStep;
                        AdaptiveStepsizeIntegrator integrator = new DormandPrince853Integrator(minStep, maxStep, tol[0], tol[1]);
                        ((AdaptiveStepsizeIntegrator)integrator).setInitialStepSize(initialStepSize);
                        propagator = new NumericalPropagator(integrator);
                        propagator.setOrbitType(orbitType);
                        propagator.setInitialState(initialState);
                        propagator.setAttitudeProvider(new LofOffset(eme2000, LOFType.VNC));
                        System.out.println("Initial orbit:");
                        System.out.println("Date:" + ((Orbit)iniOrbit).getDate());
                        System.out.println("SMA:" + ((Orbit)iniOrbit).getA() / 1000.0 + " km");
                        System.out.println("MASS:" + initialState.getMass() + " kg");
                        System.out.println("Mean:" + FastMath.toDegrees((new KeplerianOrbit(iniOrbit)).getMeanAnomaly()) + " °");
                        System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                        System.out.println("ECC:" + ((Orbit)iniOrbit).getE());
                        System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getPerigeeArgument(), Math.PI)));
                        System.out.println("Mean Motion :" + 86400.0 * initialState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                        System.out.println("Incli :" + FastMath.toDegrees(initialState.getI()));
                        System.out.println("Fram :" + initialState.getFrame());
                        System.out.println("FramP :" + propagator.getFrame());
                        System.out.println("OrbitTypeP :" + propagator.getOrbitType());
                        System.out.println("PV:" + initialState.getOrbit().getPVCoordinates());
                        System.out.println("Kep:" + initialState.getOrbit());
                        if (maneuverType.equals("Continuous")) {
                            System.out.println("Maneuver : Continuous");
                            ManeuverTriggers triggers = new DateBasedManeuverTriggers(manoeuverStartDate, durationOfManoeuver);
                            System.out.println("manoeuverStartDate :" + manoeuverStartDate);
                            ThrustPropulsionModel propulsionModel = new BasicConstantThrustPropulsionModel(THRUST, ISP, direction, "apogee-engine");
                            propagator.addForceModel(new Maneuver((AttitudeProvider)null, triggers, propulsionModel));
                        } else if (maneuverType.equals("Impulse")) {
                            System.out.println("Maneuver : Impulse");
                            EventDetector trigger = new DateDetector(manoeuverStartDate);
                            Vector3D directionImp = new Vector3D(DVx * DV, DVy * DV, DVz * DV);
                            ImpulseManeuver<EventDetector> maneuver = new ImpulseManeuver(trigger, new LofOffset(eme2000, LOFType.VNC), directionImp, ISP);
                            System.out.println(" Maneuvre :  trigger= " + trigger + " direction = " + directionImp + " ISP =" + ISP);
                            propagator.addEventDetector(maneuver);
                        } else {
                            System.out.println("Maneuver : NONE");
                        }

                        ForceModel potential = new NewtonianAttraction(MU);
                        propagator.addForceModel(potential);
                        propagator.getMultiplexer().add(stepSize, (state) -> {
                            System.out.println("Date" + state.getDate() + " e= " + state.getE() + " a= " + state.getA() + " I :" + FastMath.toDegrees(state.getI()) + " M= " + state.getMass());
                        });
                        SpacecraftState finalState = propagator.propagate(manoeuverEndDate);
                        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
                        System.out.println("Final orbit:");
                        System.out.println("Date:" + finalState.getDate());
                        System.out.println("SMA:" + finalState.getA() / 1000.0 + " km");
                        System.out.println("MeanAnomaly = " + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI)) + " deg");
                        System.out.println("manoeuverStartDate:" + manoeuverStartDate + " .");
                        System.out.println("manoeuverEndDate:" + manoeuverEndDate + " .");
                        System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                        System.out.println("ECC:" + ((Orbit)finalOrbit).getE());
                        System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI)));
                        System.out.println("Mean Motion :" + 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                        System.out.println("Incli :" + FastMath.toDegrees(finalState.getI()));
                        System.out.println("state :" + propagator.getInitialState().getDate());
                        System.out.println("PV:" + finalState.getOrbit().getPVCoordinates());
                        System.out.println("Kep:" + finalState.getOrbit());
                        System.out.println("Quaternion: at" + finalState.getDate() + " " + finalState.getAttitude().getOrientation().getRotation().getQ0() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ1() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ2() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ3());
                        FileWriter writer = new FileWriter("Result.txt", true);
                        BufferedWriter bufferedWriter = new BufferedWriter(writer);
                        bufferedWriter.newLine();
                        bufferedWriter.write("Orbital parameters post-maneuver :");
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getMass() - dryMass));
                        bufferedWriter.newLine();
                        bufferedWriter.write(((Orbit)finalOrbit).getDate().toString());
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", (finalState.getPVCoordinates().getPosition().getNorm() - 6378137.0) / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.12f", finalState.getA() / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.7f", finalState.getE()));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(finalState.getI())));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.14f", 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getKeplerianPeriod()));
                        bufferedWriter.newLine();
                        bufferedWriter.close();
                        double end = (double)System.currentTimeMillis();
                        double duration = (end - start) / 1000.0;
                        System.out.println("Execution time:" + duration);
                    } else {
                        System.out.println("Integrator : RungKutta4");
                        stepSize = 1.0;
                        ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(stepSize);
                        propagator = new NumericalPropagator(integrator);
                        propagator.setOrbitType(orbitType);
                        propagator.setInitialState(initialState);
                        propagator.setAttitudeProvider(new LofOffset(eme2000, LOFType.VNC));
                        System.out.println("Initial orbit:");
                        System.out.println("Date:" + ((Orbit)iniOrbit).getDate());
                        System.out.println("SMA:" + ((Orbit)iniOrbit).getA() / 1000.0 + " km");
                        System.out.println("MASS:" + initialState.getMass() + " kg");
                        System.out.println("Mean:" + FastMath.toDegrees((new KeplerianOrbit(iniOrbit)).getMeanAnomaly()) + " °");
                        System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                        System.out.println("ECC:" + ((Orbit)iniOrbit).getE());
                        System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)iniOrbit).getPerigeeArgument(), Math.PI)));
                        System.out.println("Mean Motion :" + 86400.0 * initialState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                        System.out.println("Incli :" + FastMath.toDegrees(initialState.getI()));
                        System.out.println("Fram :" + initialState.getFrame());
                        System.out.println("FramP :" + propagator.getFrame());
                        System.out.println("OrbitTypeP :" + propagator.getOrbitType());
                        System.out.println("PV:" + initialState.getOrbit().getPVCoordinates());
                        System.out.println("Kep:" + initialState.getOrbit());
                        if (maneuverType.equals("Continuous")) {
                            System.out.println("Maneuver : Continuous");
                            ManeuverTriggers triggers = new DateBasedManeuverTriggers(manoeuverStartDate, durationOfManoeuver);
                            System.out.println("manoeuverStartDate :" + manoeuverStartDate);
                            ThrustPropulsionModel propulsionModel = new BasicConstantThrustPropulsionModel(THRUST, ISP, direction, "apogee-engine");
                            propagator.addForceModel(new Maneuver((AttitudeProvider)null, triggers, propulsionModel));
                        } else if (maneuverType.equals("Impulse")) {
                            System.out.println("Maneuver : Impulse");
                            EventDetector trigger = new DateDetector(manoeuverStartDate);
                            Vector3D directionImp = new Vector3D(DVx * THRUST, DVy * THRUST, DVz * THRUST);
                            ImpulseManeuver<EventDetector> maneuver = new ImpulseManeuver(trigger, new LofOffset(eme2000, LOFType.VNC), directionImp, ISP);
                            System.out.println(" Maneuvre :  trigger= " + trigger + " direction = " + directionImp + " ISP =" + ISP);
                            propagator.addEventDetector(maneuver);
                        } else {
                            System.out.println("Maneuver : NONE");
                        }

                        ForceModel potential = new NewtonianAttraction(MU);
                        propagator.addForceModel(potential);
                        propagator.getMultiplexer().add(stepSize, (state) -> {
                            System.out.println("Date" + state.getDate() + " e= " + state.getE() + " a= " + state.getA() + " I :" + FastMath.toDegrees(state.getI()) + " M= " + state.getMass());
                        });
                        SpacecraftState finalState = propagator.propagate(manoeuverEndDate);
                        Orbit finalOrbit = new KeplerianOrbit(finalState.getOrbit());
                        System.out.println("Final orbit:");
                        System.out.println("Date:" + finalState.getDate());
                        System.out.println("SMA:" + finalState.getA() / 1000.0 + " km");
                        System.out.println("MeanAnomaly = " + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI)) + " deg");
                        System.out.println("manoeuverStartDate:" + manoeuverStartDate + " .");
                        System.out.println("manoeuverEndDate:" + manoeuverEndDate + " .");
                        System.out.println("RAAN:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI)));
                        System.out.println("ECC:" + ((Orbit)finalOrbit).getE());
                        System.out.println("AOP:" + FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI)));
                        System.out.println("Mean Motion :" + 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586);
                        System.out.println("Incli :" + FastMath.toDegrees(finalState.getI()));
                        System.out.println("state :" + propagator.getInitialState().getDate());
                        System.out.println("PV:" + finalState.getOrbit().getPVCoordinates());
                        System.out.println("Kep:" + finalState.getOrbit());
                        System.out.println("Quaternion: at" + finalState.getDate() + " " + finalState.getAttitude().getOrientation().getRotation().getQ0() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ1() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ2() + ", " + finalState.getAttitude().getOrientation().getRotation().getQ3());
                        FileWriter writer = new FileWriter("Result.txt", true);
                        BufferedWriter bufferedWriter = new BufferedWriter(writer);
                        bufferedWriter.newLine();
                        bufferedWriter.write("Orbital parameters post-maneuver :");
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getMass() - dryMass));
                        bufferedWriter.newLine();
                        bufferedWriter.write(((Orbit)finalOrbit).getDate().toString());
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", (finalState.getPVCoordinates().getPosition().getNorm() - 6378137.0) / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.12f", finalState.getA() / 1000.0));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.7f", finalState.getE()));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(finalState.getI())));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.4f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getRightAscensionOfAscendingNode(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getPerigeeArgument(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.8f", FastMath.toDegrees(MathUtils.normalizeAngle(((KeplerianOrbit)finalOrbit).getMeanAnomaly(), Math.PI))));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.14f", 86400.0 * finalState.getOrbit().getKeplerianMeanMotion() / 6.283185307179586));
                        bufferedWriter.newLine();
                        bufferedWriter.write(String.format("%.6f", finalState.getKeplerianPeriod()));
                        bufferedWriter.newLine();
                        bufferedWriter.close();
                        double end = (double)System.currentTimeMillis();
                        double duration = (end - start) / 1000.0;
                        System.out.println("Execution time:" + duration);
                    }
                } catch (OrekitException var92) {
                    OrekitException e = var92;
                    System.err.println(e.getLocalizedMessage());
                    System.exit(1);
                }
            } catch (IOException var93) {
                IOException e = var93;
                e.printStackTrace();
            }
        } catch (IOException var94) {
            IOException e = var94;
            e.printStackTrace();
        }

    }
}
