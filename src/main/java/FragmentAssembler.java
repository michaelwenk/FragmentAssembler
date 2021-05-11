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
        final int minMatchingSphereCount = 2;
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

        sscList = sscList.stream()
                         .filter(ssc -> {
                             final Assignment matchAssignment = Match.matchSpectra(ssc.getSpectrum(), querySpectrum, 0,
                                                                                   0, shiftTol, true, true, true);
                             System.out.println(matchAssignment.getSetAssignmentsCountWithEquivalences(0)
                                                        + " vs. "
                                                        + (ssc.getStructure()
                                                              .getAtomCount()
                                     - Utils.countHeteroAtoms(ssc.getStructure())));
                             return matchAssignment.getSetAssignmentsCountWithEquivalences(0)
                                     == ssc.getStructure()
                                           .getAtomCount()
                                     - Utils.countHeteroAtoms(ssc.getStructure());
                         })
                         .collect(Collectors.toList());

        sscList.sort((ssc1, ssc2) -> {
            final Double rmsdSSC1 = Match.calculateRMSD(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
                                                        true);
            final Double rmsdSSC2 = Match.calculateRMSD(ssc1.getSpectrum(), querySpectrum, 0, 0, shiftTol, true, true,
                                                        true);
            final int rmsdComparison = rmsdSSC1.compareTo(rmsdSSC2);

            if (rmsdComparison
                    != 0) {
                return rmsdComparison;
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

            return Integer.compare(Utils.countHeteroAtoms(ssc1.getStructure()),
                                   Utils.countHeteroAtoms(ssc2.getStructure()));
        });

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
        final Set<String> smilesSet = new HashSet<>();

        for (int i = 0; i
                < sscList.size(); i++) {
            for (int j = 0; j
                    < sscList.size(); j++) {

                //        final int i = 0;
                //        final int j = 1;
                System.out.println("\n\n--> ssc pair: "
                                           + i
                                           + " vs. "
                                           + j);
                final SSC ssc1 = sscList.get(i);
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
                        for (final SSC extendedSSC : Assembly.buildExtendedSSCList(querySpectrum, mf, maxSphere, ssc1,
                                                                                   ssc2, i, j, rootAtomIndices[0],
                                                                                   rootAtomIndices[1], shiftTol,
                                                                                   matchFactorThrs)) {
                            // do not add duplicates
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                System.out.println("SMILES: "
                                                           + smiles);
                                if (!smilesSet.contains(smiles)) {
                                    System.out.println(" --> added");
                                    extendedSSCList.add(extendedSSC);
                                    smilesSet.add(smiles);
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
                                    if (AssemblyUtils.isFinalSSC(extendedSSC, querySpectrum, shiftTol, matchFactorThrs,
                                                                 mf)) {
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
                                }
                            } catch (final CDKException e) {
                                e.printStackTrace();
                            }
                        }
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
