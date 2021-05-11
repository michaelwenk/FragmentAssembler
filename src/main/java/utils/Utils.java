package utils;

import casekit.nmr.hose.HOSECodeBuilder;
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
