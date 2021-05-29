import casekit.nmr.model.Assignment;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Match;
import model.Query;
import model.SSC;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import utils.AssemblyUtils;
import utils.Converter;
import utils.Utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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

        final int nStarts = 20;
        final int maxSphere = 4;
        final int minMatchingSphere = 1;
        final double shiftTol = 10;
        final double matchFactorThrs = 7;
        final int nThreads = 3;
        //        List<SSC> sscList = Fragmentation.buildFromNMRShiftDB("/Users/mwenk/Downloads/test.sd", new String[]{"13C"},
        //                                                              maxSphere, nThreads);
        //        List<SSC> sscList = Fragmentation.buildFromNMRShiftDB("/Users/mwenk/Downloads/nmrshiftdb2withsignals.sd",
        //                                                              new String[]{"13C"}, maxSphere, nThreads);
        List<SSC> sscList = new ArrayList<>();

        //        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList.json";
        //        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList_all_3.json";
        final String pathToJsonFile = "/Users/mwenk/Downloads/sscList_all_4.json";
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

        sscList = sscList.stream()
                         .filter(ssc -> {
                             //                             if (!ssc.getSpectrum()
                             //                                     .getSolvent()
                             //                                     .equals(querySpectrum.getSolvent())) {
                             //                                 return false;
                             //                             }
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

        final List<SSC> startSSCList = sscList.stream()
                                              .filter(ssc -> !(ssc.getUnsaturatedAtomIndices()
                                                                  .isEmpty()
                                                      && ssc.getSpectrum()
                                                            .getSignalCountWithEquivalences()
                                                      != querySpectrum.getSignalCountWithEquivalences()))
                                              .collect(Collectors.toList());
        System.out.println("\nvalid start SSC count: "
                                   + startSSCList.size());
        // prepared for possible start SSC selection
        final List<Integer> startSSCIndices = new ArrayList<>();
        //        sscList = startSSCList.subList(1, 2);
        //        startSSCIndices.add(0);
        for (int i = 0; i
                < startSSCList.size(); i++) {
            if (i
                    >= nStarts) {
                break;
            }
            startSSCIndices.add(i);
        }

        final DepictionGenerator depictionGenerator = new DepictionGenerator().withAtomColors()
                                                                              .withAtomNumbers()
                                                                              .withAromaticDisplay()
                                                                              .withSize(512, 512)
                                                                              .withFillToFit();
        try {
            for (int i = 0; i
                    < startSSCList.size()
                    && i
                    < 100; i++) {
                //                            < startSSCList.size(); i++) {
                depictionGenerator.depict(startSSCList.get(i)
                                                      .getStructure())
                                  .writeTo("/Users/mwenk/Downloads/depictions/ssc"
                                                   + i
                                                   + ".png");
            }
        } catch (final CDKException | IOException e) {
            e.printStackTrace();
        }

        //        final Map<String, SSC> solutions = Assembly.assembleWithBacktracking(sscList, querySpectrum, mf, shiftTol,
        //                                                                             matchFactorThrs, minMatchingSphere,
        //                                                                             nStarts);
        final Map<String, SSC> solutions = Assembly.assembleSequentially(sscList, querySpectrum, mf, shiftTol,
                                                                         matchFactorThrs, minMatchingSphere, nStarts);

        System.out.println("\n\nsolutions:\n"
                                   + solutions.keySet());
    }
}
