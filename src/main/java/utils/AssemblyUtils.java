package utils;

import casekit.nmr.Utils;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import model.SSC;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.*;

public class AssemblyUtils {

    public static boolean isValidExtension(final SSC ssc, final Spectrum querySpectrum, final double shiftTol,
                                           final double matchFactorThrs, final String mf,
                                           final IAtomContainer prevSubstructure) {
        return ((ssc.getStructure()
                    .getAtomCount()
                > prevSubstructure.getAtomCount())
                || (ssc.getStructure()
                       .getBondCount()
                > prevSubstructure.getBondCount()))
                && isValidSSC(ssc, querySpectrum, shiftTol, matchFactorThrs, mf);
    }

    public static boolean isValidSSC(final SSC ssc, final Spectrum querySpectrum, final double shiftTol,
                                     final double matchFactorThrs, final String mf) {

        System.out.println(isValidSubspectrum(ssc, querySpectrum, shiftTol, matchFactorThrs)
                                   + " && "
                                   + isValidSubstructure(ssc.getStructure(), mf));
        return isValidSubspectrum(ssc, querySpectrum, shiftTol, matchFactorThrs)
                && isValidSubstructure(ssc.getStructure(), mf);
    }

    public static boolean isValidSubspectrum(final SSC ssc, final Spectrum querySpectrum, final double shiftTol,
                                             final double matchFactorThrs) {
        if ((ssc.getSpectrum()
                == null)
                || (ssc.getSpectrum()
                       .getSignalCount()
                > querySpectrum.getSignalCount())) {
            System.out.println("-> subspectrum == null or signal count subspectrum > signal count query spectrum!!!");
            return false;
        }
        final Assignment matchAssignments = Match.matchSpectra(ssc.getSpectrum(), querySpectrum, 0, 0, shiftTol, true,
                                                               true, true);
        // filter for unset assignments
        if (matchAssignments.getSetAssignmentsCount(0)
                == 0) {
            return false;
        }
        // @TODO add/enable filter for intensities
        final Double avgDev = Match.calculateAverageDeviation(ssc.getSpectrum(), querySpectrum, 0, 0, shiftTol, true,
                                                              true, true);
        if (avgDev
                == null
                || avgDev
                > matchFactorThrs) {
            System.out.println("-> match factor not allowed!!!");
            return false;
        }

        return true;
    }

    public static boolean isValidSubstructure(final IAtomContainer substructure, final String mf) {
        double bondOrderSum;
        IAtom atom;
        for (int i = 0; i
                < substructure.getAtomCount(); i++) {
            bondOrderSum = Utils.getBondOrderSum(substructure, i, true);
            atom = substructure.getAtom(i);
            // @TODO include different valencies: N3, N5, S2, S4, S6 etc.
            // -1 for cases with heterocyclic aromatics, like the N in the small aromatic ring in coffein if we want to add the bond to the CH3 group
            if (atom.isAromatic()
                    && (atom.getSymbol()
                            .equals("N")
                    || atom.getSymbol()
                           .equals("S")
                    || atom.getSymbol()
                           .equals("P"))) {
                //            System.out.print("[ -1 ]");
                bondOrderSum -= 1;
            }
            if (bondOrderSum
                    > substructure.getAtom(i)
                                  .getValency()) {
                return false;
            }
        }
        return compareWithMolecularFormula(substructure, mf);
    }

    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum, final double shiftTol,
                                     final double matchFactorThrs, final String mf) {
        System.out.println("\nunsaturated atoms left? -> "
                                   + !ssc.getUnsaturatedAtomIndices()
                                         .isEmpty()
                                   + " -> "
                                   + ssc.getUnsaturatedAtomIndices());
        if (!ssc.getUnsaturatedAtomIndices()
                .isEmpty()) {
            return false;
        }
        System.out.println("query spectrum size reached? -> "
                                   + ssc.getSpectrum()
                                        .getSignalCountWithEquivalences()
                                   + " == "
                                   + querySpectrum.getSignalCountWithEquivalences()
                                   + " -> "
                                   + (ssc.getSpectrum()
                                         .getSignalCountWithEquivalences()
                == querySpectrum.getSignalCountWithEquivalences()));
        if ((ssc.getSpectrum()
                .getSignalCountWithEquivalences()
                != querySpectrum.getSignalCountWithEquivalences())) {
            return false;
        }
        System.out.println("isValidSpectrum? -> "
                                   + AssemblyUtils.isValidSSC(ssc, querySpectrum, shiftTol, matchFactorThrs, mf));
        if (!isValidSSC(ssc, querySpectrum, shiftTol, matchFactorThrs, mf)) {
            return false;
        }
        final IMolecularFormula molecularFormula = casekit.nmr.utils.Utils.getMolecularFormulaFromString(mf);
        System.out.println("MF vs. atom count: "
                                   + ((molecularFormula
                != null)
                && (MolecularFormulaManipulator.getAtomCount(molecularFormula)
                != MolecularFormulaManipulator.getAtomCount(
                MolecularFormulaManipulator.getMolecularFormula(ssc.getStructure())))));
        if ((molecularFormula
                != null)
                && (MolecularFormulaManipulator.getAtomCount(molecularFormula)
                != MolecularFormulaManipulator.getAtomCount(
                MolecularFormulaManipulator.getMolecularFormula(ssc.getStructure())))) {
            return false;
        }
        try {
            Kekulization.kekulize(ssc.getStructure());
            System.out.println("kekulization? -> true");
        } catch (final CDKException e) {
            System.out.println("kekulization? -> false");
            return false;
        }

        return true;
    }

    public static int getMaximumMatchingSphereHOSECode(final SSC ssc1, final SSC ssc2, final int atomIndexSSC1,
                                                       final int atomIndexSSC2) {
        int maxMatchingSphere = -1;

        if (!Utils.checkIndexInAtomContainer(ssc1.getStructure(), atomIndexSSC1)
                || !Utils.checkIndexInAtomContainer(ssc2.getStructure(), atomIndexSSC2)) {
            return maxMatchingSphere;
        }

        final List<String> HOSECodeSpheresSSC1 = casekit.nmr.hose.Utils.splitHOSECodeIntoSpheres(ssc1.getHoseCodes()
                                                                                                     .get(atomIndexSSC1));
        final List<String> HOSECodeSpheresSSC2 = casekit.nmr.hose.Utils.splitHOSECodeIntoSpheres(ssc2.getHoseCodes()
                                                                                                     .get(atomIndexSSC2));
        for (int s = 0; s
                < Integer.min(HOSECodeSpheresSSC1.size(), HOSECodeSpheresSSC2.size()); s++) {
            if (!HOSECodeSpheresSSC1.get(s)
                                    .equals(HOSECodeSpheresSSC2.get(s))) {
                break;
            }
            maxMatchingSphere = s;
        }

        return maxMatchingSphere;
    }

    public static Map<Integer, List<Integer[]>> getOverlaps(final SSC ssc1, final SSC ssc2,
                                                            final int minMatchingSphereCount) {
        final Map<Integer, List<Integer[]>> overlapsInSpheres = new HashMap<>();
        int maxMatchingSphere;
        for (int i = 0; i
                < ssc1.getStructure()
                      .getAtomCount(); i++) {
            for (int j = 0; j
                    < ssc2.getStructure()
                          .getAtomCount(); j++) {
                maxMatchingSphere = getMaximumMatchingSphereHOSECode(ssc1, ssc2, i, j);
                if (maxMatchingSphere
                        < minMatchingSphereCount) {
                    continue;
                }
                // create new key for the found max. matching sphere if it's not existing
                overlapsInSpheres.putIfAbsent(maxMatchingSphere, new ArrayList<>());
                // insert matching atom pair into mappings hashmap with found max. matching sphere
                overlapsInSpheres.get(maxMatchingSphere)
                                 .add(new Integer[]{i, j});
            }
        }

        return overlapsInSpheres;
    }

    public static boolean compareWithMolecularFormula(final IAtomContainer structure, final String mf) {
        if (mf
                != null
                && !mf.trim()
                      .isEmpty()) {
            for (final String atomType : Utils.getAtomTypesInAtomContainer(structure)) {
                if (Utils.getAtomTypeIndicesByElement(structure, atomType)
                         .size()
                        > MolecularFormulaManipulator.getElementCount(
                        casekit.nmr.utils.Utils.getMolecularFormulaFromString(mf), atomType)) {
                    return false;
                }
            }
        }

        return true;
    }

    public static void addSignalToSSC(final Spectrum spectrum, final Assignment assignment, final Signal signalToAdd,
                                      final int atomIndex) {
        if (signalToAdd
                == null) {
            return;
        }

        // just to be sure that we take the right signal if equivalences are present
        final List<Integer> closestSignalList = spectrum.pickByClosestShift(signalToAdd.getShift(0), 0, 0.0);
        closestSignalList.retainAll(spectrum.pickByMultiplicity(signalToAdd.getMultiplicity()));
        if (closestSignalList.isEmpty()) {
            spectrum.addSignal(signalToAdd);
            assignment.addAssignment(0, new int[]{atomIndex});
        } else {
            final int signalIndex = closestSignalList.get(0);
            if (Arrays.stream(assignment.getAssignment(0, signalIndex))
                      .noneMatch(equiv -> equiv
                              == atomIndex)) {
                assignment.addAssignmentEquivalence(0, signalIndex, atomIndex);
            }
        }
    }
}
