package utils;

import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class Utils {

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

    public static List<String> buildHOSECodes(final IAtomContainer ac, final int maxSphere) {
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

    public static void sortSSCList(final List<SSC> sscList, final Spectrum querySpectrum, final double shiftTol) {
        sscList.sort((ssc1, ssc2) -> {
            final Double rmsdSSC1 = Match.calculateRMSD(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
                                                        true);
            final Double rmsdSSC2 = Match.calculateRMSD(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
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

            // @TODO consider number of rings ?
            // @TODO consider number of open sites ?

            return Integer.compare(Utils.countHeteroAtoms(ssc1.getStructure()),
                                   Utils.countHeteroAtoms(ssc2.getStructure()));
        });
    }

    public static int countHeteroAtoms(final IAtomContainer structure) {
        int heteroAtomCount = 0;
        for (int i = 0; i
                < structure.getAtomCount(); i++) {
            if (
                //                    !structure.getAtom(i)
                //                          .getSymbol()
                //                          .equals("H")
                //                    &&
                    !structure.getAtom(i)
                              .getSymbol()
                              .equals("C")) {
                heteroAtomCount++;
            }
        }
        return heteroAtomCount;
    }
}
