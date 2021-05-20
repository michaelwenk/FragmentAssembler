import casekit.nmr.dbservice.NMRShiftDB;
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
        //        List<SSC> sscList = buildFromNMRShiftDB("/Users/mwenk/Downloads/test.sd", new String[]{"13C"}, maxSphere, 2);
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

        final int i = 0;
        SSC ssc1 = sscList.get(i);
        for (int j = 0; j
                < sscList.size(); j++) {
            if (i
                    == j) {
                continue;
            }
            System.out.println("\n\n--> ssc pair: "
                                       + i
                                       + " vs. "
                                       + j);
            final SSC ssc2 = sscList.get(j);
            final Map<Integer, List<Integer[]>> overlaps = AssemblyUtils.getOverlaps(ssc1, ssc2,
                                                                                     minMatchingSphereCount);
            final List<SSC> extendedSSCList = new ArrayList<>();
            String smiles;
            for (final Map.Entry<Integer, List<Integer[]>> integerListEntry : overlaps.entrySet()) {
                System.out.println("---> sphere match: "
                                           + integerListEntry.getKey());
                for (final Integer[] rootAtomIndices : integerListEntry.getValue()) {
                    System.out.println("\n--> root atom pair: "
                                               + Arrays.toString(rootAtomIndices));
                    for (final SSC extendedSSC : Assembly.assemblyCore(querySpectrum, mf, maxSphere, ssc1, ssc2,
                                                                       rootAtomIndices[0], rootAtomIndices[1], shiftTol,
                                                                       matchFactorThrs)) {
                        // do not add duplicates
                        try {
                            smiles = smilesGenerator.create(extendedSSC.getStructure());
                            System.out.println("SMILES: "
                                                       + smiles);
                            System.out.println(" --> added");
                            extendedSSCList.add(extendedSSC);
                            try {
                                depictionGenerator.depict(extendedSSC.getStructure())
                                                  .writeTo("/Users/mwenk/Downloads/depictions/extended_"
                                                                   + i
                                                                   + "-"
                                                                   + j
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
                                                                       + i
                                                                       + "-"
                                                                       + j
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
                System.out.println("SIZE: "
                                           + extendedSSCList.size());
                Utils.sortSSCList(extendedSSCList, querySpectrum, shiftTol);
                ssc1 = extendedSSCList.get(0)
                                      .buildClone();
            }
        }
    }

    public static List<SSC> buildFromNMRShiftDB(final String pathToNMRShiftDB, final String[] nuclei,
                                                final int maxSphere, final int nThreads) {
        final List<SSC> sscList = new ArrayList<>();
        List<SSC> sscListTemp = new ArrayList<>();
        for (int m = 2; m
                <= maxSphere; m++) {
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
        }
        System.out.println("Building SSC done!!!");

        return sscList;
    }
}
