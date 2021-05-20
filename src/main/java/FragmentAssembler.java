import casekit.nmr.model.Assignment;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import model.Query;
import model.SSC;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import utils.AssemblyUtils;
import utils.Converter;
import utils.Utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class FragmentAssembler {

    private DepictionGenerator depictionGenerator;

    public static void main(final String[] args) {
        final List<Query> queries = Converter.parseQueriesFile("/Users/mwenk/Downloads/queries.txt");
        System.out.println("#queries: "
                                   + queries.size());

        final Spectrum querySpectrum = queries.get(0)
                                              .getQuerySpectrum();
        final String mf = queries.get(0)
                                 .getMf();
        System.out.println(querySpectrum);
        System.out.println(mf);

        final int maxSphere = 3;
        final int minMatchingSphereCount = 1;
        final double shiftTol = 2;
        final double matchFactorThrs = 2;
        //        List<SSC> sscList = Fragmentation.buildFromNMRShiftDB("/Users/mwenk/Downloads/test.sd", new String[]{"13C"},
        //                                                              maxSphere, 2);
        List<SSC> sscList = new ArrayList<>();

        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList.json";
        //        Converter.sscListToJSONFile(sscList, pathToJsonFile);

        try {
            sscList = Converter.jsonFileToSscList(pathToJsonFile, 2);
        } catch (final InterruptedException e) {
            e.printStackTrace();
        } catch (final FileNotFoundException e) {
            sscList = new ArrayList<>();
            e.printStackTrace();
        }

        sscList = sscList.parallelStream()
                         .filter(ssc -> {
                             final Assignment matchAssignment = Match.matchSpectra(ssc.getSpectrum(), querySpectrum, 0,
                                                                                   0, shiftTol, true, true, true);
                             return matchAssignment.getSetAssignmentsCountWithEquivalences(0)
                                     == ssc.getStructure()
                                           .getAtomCount()
                                     - Utils.countHeteroAtoms(ssc.getStructure());
                         })
                         .collect(Collectors.toList());

        Utils.sortSSCList(sscList, querySpectrum, shiftTol);

        System.out.println("\ntotal SSC count: "
                                   + sscList.size());

        final DepictionGenerator depictionGenerator = new DepictionGenerator().withAtomColors()
                                                                              .withAtomNumbers()
                                                                              .withAromaticDisplay()
                                                                              .withSize(512, 512)
                                                                              .withFillToFit();
        //        try {
        //            for (int i = 0; i
        //                    < sscList.size(); i++) {
        //                depictionGenerator.depict(sscList.get(i)
        //                                                 .getStructure())
        //                                  .writeTo("/Users/mwenk/Downloads/depictions/ssc"
        //                                                   + i
        //                                                   + ".png");
        //            }
        //        } catch (final CDKException | IOException e) {
        //            e.printStackTrace();
        //        }

        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        final Map<String, SSC> solutions = new HashMap<>();
        final Stack<Object[]> intermediates = new Stack<>();

        final int startSSCIndex = 0;
        intermediates.push(new Object[]{sscList.get(startSSCIndex).buildClone(), startSSCIndex});

        Object[] stackItem;
        SSC intermediate;
        int currentSSCIndex;
        SSC currentSSC;
        List<SSC> extendedSSCList;
        String smiles;
        Map<Integer, List<Integer[]>> overlaps;

        while (!intermediates.empty()) {
            stackItem = intermediates.pop();
            intermediate = (SSC) stackItem[0];
            currentSSCIndex = ((int) stackItem[1])
                    + 1;
            if (currentSSCIndex
                    >= sscList.size()) {
                continue;
            }

            //            if (currentSSCIndex
            //                    == startSSCIndex) {
            //                continue;
            //            }
            System.out.println("\n\n--> ssc pair: "
                                       + startSSCIndex
                                       + " vs. "
                                       + currentSSCIndex);
            currentSSC = sscList.get(currentSSCIndex);
            overlaps = AssemblyUtils.getOverlaps(intermediate, currentSSC, minMatchingSphereCount);
            extendedSSCList = new ArrayList<>();
            for (final Map.Entry<Integer, List<Integer[]>> integerListEntry : overlaps.entrySet()) {
                System.out.println("---> sphere match: "
                                           + integerListEntry.getKey());
                for (final Integer[] rootAtomIndices : integerListEntry.getValue()) {
                    System.out.println("\n--> root atom pair: "
                                               + Arrays.toString(rootAtomIndices));
                    for (final SSC extendedSSC : Assembly.assemblyCore(querySpectrum, mf, maxSphere,
                                                                       intermediate.buildClone(), currentSSC,
                                                                       rootAtomIndices[0], rootAtomIndices[1], shiftTol,
                                                                       matchFactorThrs)) {
                        extendedSSCList.add(extendedSSC);
                        try {
                            smiles = smilesGenerator.create(extendedSSC.getStructure());
                            try {
                                depictionGenerator.depict(extendedSSC.getStructure())
                                                  .writeTo("/Users/mwenk/Downloads/depictions/extended_"
                                                                   + startSSCIndex
                                                                   + "-"
                                                                   + currentSSCIndex
                                                                   + "_"
                                                                   + rootAtomIndices[0]
                                                                   + "-"
                                                                   + rootAtomIndices[1]
                                                                   + ".png");
                            } catch (final IOException | CDKException e) {
                                e.printStackTrace();
                            }
                            if (AssemblyUtils.isFinalSSC(extendedSSC, querySpectrum, shiftTol, matchFactorThrs, mf)) {
                                solutions.putIfAbsent(smiles, extendedSSC);
                                try {
                                    depictionGenerator.depict(extendedSSC.getStructure())
                                                      .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                                       + startSSCIndex
                                                                       + "-"
                                                                       + currentSSCIndex
                                                                       + "_"
                                                                       + rootAtomIndices[0]
                                                                       + "-"
                                                                       + rootAtomIndices[1]
                                                                       + ".png");
                                } catch (final IOException | CDKException e) {
                                    e.printStackTrace();
                                }
                            }
                        } catch (final CDKException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
            if (!extendedSSCList.isEmpty()) {
                Utils.sortSSCList(extendedSSCList, querySpectrum, shiftTol);
                Collections.reverse(extendedSSCList);
                for (final SSC newIntermediate : extendedSSCList) {
                    intermediates.push(new Object[]{newIntermediate.buildClone(), currentSSCIndex});
                }
            } else {
                intermediates.push(new Object[]{intermediate, currentSSCIndex});
            }
        }

        System.out.println("\n\nsolutions:\n"
                                   + solutions.keySet());
    }
}
