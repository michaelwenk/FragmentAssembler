package utils;

import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.*;

public class Utils {

    final private static SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);

    public static List<Integer> getUnsaturatedAtomIndices(final IAtomContainer ac) {
        final List<Integer> unsaturatedAtomIndices = new ArrayList<>();
        for (int i = 0; i
                < ac.getAtomCount(); i++) {
            // set the indices of unsaturated atoms in substructure
            if (!casekit.nmr.Utils.isSaturated(ac, i)) {
                unsaturatedAtomIndices.add(i);
            }
        }
        return unsaturatedAtomIndices;
    }

    public static List<String> buildHOSECodes(final IAtomContainer ac, final Integer maxSphere) {
        final List<String> hoseCodes = new ArrayList<>();
        for (int i = 0; i
                < ac.getAtomCount(); i++) {
            try {
                hoseCodes.add(HOSECodeBuilder.buildHOSECode(ac, i, maxSphere, false));
            } catch (final CDKException e) {
                e.printStackTrace();
                return null;
            }
        }
        return hoseCodes;
    }

    public static Integer findKeyInMap(final Map<Integer, Integer> map, final int value) {
        final Optional<Integer> keyOptional = map.entrySet()
                                                 .stream()
                                                 .filter(entry -> entry.getValue()
                                                         == value)
                                                 .map(Map.Entry::getKey)
                                                 .findFirst();
        return keyOptional.orElse(null);
    }

    public static void removeDuplicatesFromSSCList(final List<SSC> sscList, final Set<String> smilesSet) {
        final List<SSC> toRemove = new ArrayList<>();
        String smiles;
        for (final SSC ssc : sscList) {
            try {
                smiles = smilesGenerator.create(ssc.getStructure());
                if (!smilesSet.contains(smilesSet)) {
                    smilesSet.add(smiles);
                } else {
                    toRemove.add(ssc);
                }
            } catch (final CDKException e) {
                e.printStackTrace();
            }
        }
        sscList.removeAll(toRemove);
    }

    public static void sortSSCList(final List<SSC> sscList, final Spectrum querySpectrum, final double shiftTol) {
        sscList.sort((ssc1, ssc2) -> {

            final Assignment matchAssignmentSSC1 = Match.matchSpectra(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol,
                                                                      true, true, true);
            final Assignment matchAssignmentSSC2 = Match.matchSpectra(ssc2.getSpectrum(), querySpectrum, 0, 0, shiftTol,
                                                                      true, true, true);
            final int matchedSignalCountComparison = -1
                    * Integer.compare(matchAssignmentSSC1.getSetAssignmentsCountWithEquivalences(0),
                                      matchAssignmentSSC2.getSetAssignmentsCountWithEquivalences(0));
            if (matchedSignalCountComparison
                    != 0) {
                return matchedSignalCountComparison;
            }

            final int atomCountComparison = -1
                    * Integer.compare(ssc1.getStructure()
                                          .getAtomCount(), ssc2.getStructure()
                                                               .getAtomCount());

            final Double rmsdSSC1 = Match.calculateRMSD(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
                                                        true);
            final Double rmsdSSC2 = Match.calculateRMSD(ssc2.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
                                                        true);
            if (rmsdSSC1
                    != null
                    && rmsdSSC2
                    != null) {
                final int rmsdComparison = rmsdSSC1.compareTo(rmsdSSC2);

                if (rmsdComparison
                        != 0) {
                    return rmsdComparison;
                }
            }

            if (atomCountComparison
                    != 0) {
                return atomCountComparison;
            }
            // consider number of open sites
            final int openSitesCountComparison = Integer.compare(ssc1.getUnsaturatedAtomIndices()
                                                                     .size(), ssc2.getUnsaturatedAtomIndices()
                                                                                  .size());
            if (openSitesCountComparison
                    != 0) {
                return openSitesCountComparison;
            }

            return Integer.compare(Utils.countHeteroAtoms(ssc1.getStructure()),
                                   Utils.countHeteroAtoms(ssc2.getStructure()));
        });
    }

    public static int countHeteroAtoms(final IAtomContainer structure) {
        int heteroAtomCount = 0;
        for (int i = 0; i
                < structure.getAtomCount(); i++) {
            if (!structure.getAtom(i)
                          .getSymbol()
                          .equals("C")) {
                heteroAtomCount++;
            }
        }
        return heteroAtomCount;
    }

    public static Map<String, Double[]> buildHOSECodeShiftStatistics(final List<SSC> sscList, final Integer maxSphere) {
        final Map<String, List<Double>> hoseCodeShifts = new HashMap<>();
        Signal signal;
        String hoseCode;
        for (final SSC ssc : sscList) {
            for (int i = 0; i
                    < ssc.getStructure()
                         .getAtomCount(); i++) {
                signal = ssc.getSpectrum()
                            .getSignal(ssc.getAssignment()
                                          .getIndex(0, i));
                if (signal
                        != null) {
                    try {
                        if (maxSphere
                                != null) {
                            for (int sphere = 0; sphere
                                    <= maxSphere; sphere++) {
                                hoseCode = HOSECodeBuilder.buildHOSECode(ssc.getStructure(), i, sphere, false);
                                hoseCodeShifts.putIfAbsent(hoseCode, new ArrayList<>());
                                hoseCodeShifts.get(hoseCode)
                                              .add(signal.getShift(0));
                            }
                        }
                        hoseCode = HOSECodeBuilder.buildHOSECode(ssc.getStructure(), i, null, false);
                        hoseCodeShifts.putIfAbsent(hoseCode, new ArrayList<>());
                        hoseCodeShifts.get(hoseCode)
                                      .add(signal.getShift(0));
                    } catch (final CDKException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        final Map<String, Double[]> hoseCodeShiftStatistics = new HashMap<>();
        for (final Map.Entry<String, List<Double>> shifts : hoseCodeShifts.entrySet()) {
            hoseCodeShiftStatistics.put(shifts.getKey(), new Double[]{Double.valueOf(shifts.getValue()
                                                                                           .size()),
                                                                      Collections.min(shifts.getValue()),
                                                                      casekit.nmr.Utils.getRMS(shifts.getValue()),
                                                                      casekit.nmr.Utils.getMedian(shifts.getValue()),
                                                                      Collections.max(shifts.getValue())});
        }
        return hoseCodeShiftStatistics;
    }
}
