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

        final int nStarts = 1;
        final int maxSphere = 3;
        final int minMatchingSphereCount = 1;
        final double shiftTol = 5;
        final double matchFactorThrs = 4;
        final int nThreads = 3;
        //        List<SSC> sscList = Fragmentation.buildFromNMRShiftDB("/Users/mwenk/Downloads/test.sd", new String[]{"13C"},
        //                                                              maxSphere, nThreads);
        //        List<SSC> sscList = Fragmentation.buildFromNMRShiftDB("/Users/mwenk/Downloads/nmrshiftdb2withsignals.sd",
        //                                                              new String[]{"13C"}, maxSphere, nThreads);
        List<SSC> sscList = new ArrayList<>();

        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList.json";
        //        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList_all_3.json";
        //        Converter.sscListToJSONFile(sscList, pathToJsonFile);

        try {
            sscList = Converter.jsonFileToSscList(pathToJsonFile, nThreads);
            System.out.println("\ntotal SSC count: "
                                       + sscList.size());
        } catch (final InterruptedException | FileNotFoundException e) {
            sscList = new ArrayList<>();
            e.printStackTrace();
        }

        //        Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        //        hoseCodeShiftStatistics = Utils.buildHOSECodeShiftStatistics(sscList, maxSphere);
        //        final String pathToHOSECodeShiftStatistics = "/Users/mwenk/Downloads/jsonToHOSECodeShiftStatistics_3.json";
        //        Converter.hoseCodeShiftStatisticsToJSONFile(hoseCodeShiftStatistics, pathToHOSECodeShiftStatistics);
        //        try {
        //            hoseCodeShiftStatistics = Converter.jsonFileToHOSECodeShiftStatistics(pathToHOSECodeShiftStatistics);
        //        } catch (final FileNotFoundException e) {
        //            hoseCodeShiftStatistics = new HashMap<>();
        //            e.printStackTrace();
        //        }

        sscList = sscList.parallelStream()
                         .filter(ssc -> {
                             if (!ssc.getSpectrum()
                                     .getSolvent()
                                     .equals(querySpectrum.getSolvent())) {
                                 return false;
                             }
                             if (!AssemblyUtils.compareWithMolecularFormula(ssc.getStructure(), mf)) {
                                 return false;
                             }

                             if (ssc.getSpectrum()
                                    .getSignalCountWithEquivalences()
                                     > querySpectrum.getSignalCountWithEquivalences()) {
                                 return false;
                             }

                             final Double averageDeviation = Match.calculateAverageDeviation(ssc.getSpectrum(),
                                                                                             querySpectrum, 0, 0,
                                                                                             shiftTol, true, true,
                                                                                             true);
                             if (averageDeviation
                                     != null
                                     && averageDeviation
                                     > matchFactorThrs) {
                                 return false;
                             }

                             final Assignment matchAssignment = Match.matchSpectra(ssc.getSpectrum(), querySpectrum, 0,
                                                                                   0, shiftTol, true, true, true);
                             return matchAssignment.getSetAssignmentsCountWithEquivalences(0)
                                     == ssc.getStructure()
                                           .getAtomCount()
                                     - Utils.countHeteroAtoms(ssc.getStructure());
                         })
                         .collect(Collectors.toList());
        Utils.sortSSCList(sscList, querySpectrum, shiftTol);

        System.out.println("\nfiltered SSC count: "
                                   + sscList.size());

        final DepictionGenerator depictionGenerator = new DepictionGenerator().withAtomColors()
                                                                              .withAtomNumbers()
                                                                              .withAromaticDisplay()
                                                                              .withSize(512, 512)
                                                                              .withFillToFit();
        //        try {
        //            for (int i = 0; i
        //                    < sscList.size()
        //                    && i
        //                    < 100; i++) {
        //                //                            < sscList.size(); i++) {
        //                depictionGenerator.depict(sscList.get(i)
        //                                                 .getStructure())
        //                                  .writeTo("/Users/mwenk/Downloads/depictions/ssc"
        //                                                   + i
        //                                                   + ".png");
        //            }
        //        } catch (final CDKException | IOException e) {
        //            e.printStackTrace();
        //        }

        // prepared for possible start SSC selection
        final List<Integer> startSSCIndices = new ArrayList<>();
        for (int i = 0; i
                < nStarts; i++) {
            startSSCIndices.add(i);
        }

        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        final Map<String, SSC> solutions = new HashMap<>();
        Stack<Object[]> intermediates;
        Object[] stackItem;
        SSC intermediate;
        int currentSSCIndex, startSSCIndex;
        List<SSC> extendedSSCList;
        String smiles;
        boolean stop;
        for (final int i : startSSCIndices) {
            intermediates = new Stack<>();
            intermediates.push(new Object[]{sscList.get(i).buildClone(), i, 0});

            while (!intermediates.empty()) {
                stackItem = intermediates.pop();
                intermediate = (SSC) stackItem[0];
                startSSCIndex = (int) stackItem[1];
                currentSSCIndex = ((int) stackItem[2]);

                if (currentSSCIndex
                        >= sscList.size()) {
                    continue;
                }
                stop = false;
                extendedSSCList = new ArrayList<>();
                for (int j = currentSSCIndex; j
                        < sscList.size(); j++) {
                    if (j
                            == startSSCIndex) {
                        continue;
                    }
                    //                    System.out.println("\n\n--> ssc pair: "
                    //                                               + startSSCIndex
                    //                                               + " vs. "
                    //                                               + j);
                    extendedSSCList.addAll(
                            Assembly.assemblyOverlaps(intermediate, sscList.get(j), querySpectrum, mf, null,//maxSphere,
                                                      shiftTol, matchFactorThrs, minMatchingSphereCount));
                }
                if (!extendedSSCList.isEmpty()) {
                    Utils.sortSSCList(extendedSSCList, querySpectrum, shiftTol);
                    Utils.removeDuplicatesFromSSCList(sscList, new HashSet<>());
                    Collections.reverse(extendedSSCList);
                    for (final SSC extendedSSC : extendedSSCList) {
                        if (AssemblyUtils.isFinalSSC(extendedSSC, querySpectrum, shiftTol, matchFactorThrs, mf)) {
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                solutions.putIfAbsent(smiles, extendedSSC);
                                try {
                                    depictionGenerator.depict(extendedSSC.getStructure())
                                                      .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                                       + startSSCIndex
                                                                       + "-"
                                                                       + currentSSCIndex
                                                                       + "_"
                                                                       + smiles
                                                                       + ".png");
                                } catch (final IOException | CDKException e) {
                                    e.printStackTrace();
                                }
                                stop = true;
                            } catch (final CDKException e) {
                                e.printStackTrace();
                            }
                        } else {
                            //                            try {
                            //                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                            //                                depictionGenerator.depict(extendedSSC.getStructure())
                            //                                                  .writeTo("/Users/mwenk/Downloads/depictions/extended_"
                            //                                                                   + startSSCIndex
                            //                                                                   + "-"
                            //                                                                   + currentSSCIndex
                            //                                                                   + "_"
                            //                                                                   + smiles
                            //                                                                   + "_"
                            //                                                                   + extendedSSC.getUnsaturatedAtomIndices()
                            //                                                                   + ".png");
                            //                            } catch (final IOException | CDKException e) {
                            //                                e.printStackTrace();
                            //                            }

                            intermediates.push(new Object[]{extendedSSC.buildClone(), startSSCIndex, currentSSCIndex
                                    + 1});
                        }
                    }
                } else {
                    intermediates.push(new Object[]{intermediate, startSSCIndex, currentSSCIndex
                            + 1});
                }
                if (stop) {
                    break;
                }
            }
        }

        System.out.println("\n\nsolutions:\n"
                                   + solutions.keySet());
    }
}
