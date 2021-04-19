import casekit.nmr.dbservice.NMRShiftDB;
import casekit.nmr.model.Spectrum;
import model.Query;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import utils.AssemblyUtils;
import utils.Converter;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class FragmentAssembler {

    public static void main(final String[] args) {
        final List<Query> queries = Converter.parseQueriesFile("/Users/mwenk/Downloads/queries.txt");
        System.out.println(queries.size());

        final Spectrum querySpectrum = queries.get(0)
                                              .getQuerySpectrum();
        final String mf = queries.get(0)
                                 .getMf();
        System.out.println(querySpectrum);
        System.out.println(mf);

        final int maxSphere = 3;
        final double shiftTol = 2;
        final double matchFactorThrs = 2;
        List<SSC> sscList = new ArrayList<>(); //buildFromNMRShiftDB("/Users/mwenk/Downloads/test.sd", new String[]{"13C"}, maxSphere, 2);

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

        System.out.println("\ntotal SSC count: "
                                   + sscList.size());

        //        final DepictionGenerator depictionGenerator = new DepictionGenerator().withAtomColors()
        //                                                                              .withAtomNumbers()
        //                                                                              .withAromaticDisplay()
        //                                                                              .withSize(512, 512)
        //                                                                              .withFillToFit();
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


        for (int i = 0; i
                < sscList.size(); i++) {
            for (int j = 0; j
                    < sscList.size(); j++) {

                //        final int i = 15;
                //        final int j = 20;
                System.out.println("\n\n--> ssc pair: "
                                           + i
                                           + " vs. "
                                           + j);
                final SSC ssc1 = sscList.get(i);
                final SSC ssc2 = sscList.get(j);
                System.out.println(ssc1.getSpectrum());
                System.out.println(ssc2.getSpectrum());
                final Map<Integer, List<Integer[]>> overlaps = AssemblyUtils.getOverlaps(ssc1, ssc2, 2);
                for (final Map.Entry<Integer, List<Integer[]>> integerListEntry : overlaps.entrySet()) {
                    System.out.println("---> sphere match: "
                                               + integerListEntry.getKey());
                    for (final Integer[] rootAtomIndices : integerListEntry.getValue()) {
                        System.out.println("\n--> root atom pair: "
                                                   + Arrays.toString(rootAtomIndices));
                        Assembly.buildFromSSCs(querySpectrum, mf, maxSphere, ssc1, ssc2, i, j, rootAtomIndices[0],
                                               rootAtomIndices[1], shiftTol, matchFactorThrs);
                    }
                }
            }
        }

    }


    public static List<SSC> buildFromNMRShiftDB(final String pathToNMRShiftDB, final String[] nuclei,
                                                final int maxSphere, final int nThreads) {
        final List<SSC> sscList = new ArrayList<>();
        List<SSC> sscListTemp = new ArrayList<>();
        //        for (int m = 2; m
        //                <= maxSphere; m++) {
        final int m = maxSphere;
        System.out.println("Building SSC for "
                                   + m
                                   + "-spheres...");

        try {
            sscListTemp = Fragmentation.buildSSCCollection(
                    NMRShiftDB.getDataSetsFromNMRShiftDB(pathToNMRShiftDB, nuclei), m, nThreads);
        } catch (final InterruptedException | FileNotFoundException | CDKException e) {
            e.printStackTrace();
        }

        System.out.println("SSC for "
                                   + m
                                   + "-spheres built!!!");
        System.out.println("-> #SSC in SSC library: "
                                   + sscListTemp.size());
        sscList.addAll(sscListTemp);
        //        }
        System.out.println("Building SSC done!!!");

        return sscList;
    }
}
