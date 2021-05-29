import casekit.nmr.Utils;
import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.hose.model.ConnectionTree;
import casekit.nmr.hose.model.ConnectionTreeNode;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import model.SSC;
import org.openscience.cdk.Bond;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import utils.AssemblyUtils;

import java.io.IOException;
import java.util.*;

public class Assembly {

    private static final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
    private static final DepictionGenerator depictionGenerator = new DepictionGenerator().withAtomColors()
                                                                                         .withAtomNumbers()
                                                                                         .withAromaticDisplay()
                                                                                         .withSize(512, 512)
                                                                                         .withFillToFit();

    public static List<SSC> assemblyOverlaps(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum,
                                             final String mf, final Integer maxSphere, final double shiftTol,
                                             final double matchFactorThrs, final int minMatchingSphere) {
        final Map<Integer, List<Integer[]>> overlaps = AssemblyUtils.getOverlaps(ssc1, ssc2, minMatchingSphere);
        final List<SSC> extendedSSCList = new ArrayList<>();
        for (final Map.Entry<Integer, List<Integer[]>> integerListEntry : overlaps.entrySet()) {
            //            System.out.println("---> sphere match: "
            //                                       + integerListEntry.getKey());
            for (final Integer[] rootAtomIndices : integerListEntry.getValue()) {
                //                System.out.println("\n--> root atom pair: "
                //                                           + Arrays.toString(rootAtomIndices));
                extendedSSCList.addAll(
                        Assembly.assemblyCore(querySpectrum, mf, maxSphere, ssc1, ssc2, rootAtomIndices[0],
                                              rootAtomIndices[1], shiftTol, matchFactorThrs));
            }
        }

        return extendedSSCList;
    }

    public static List<SSC> assemblyCore(final Spectrum querySpectrum, final String mf, final Integer maxSphere,
                                         final SSC ssc1, final SSC ssc2, final int rootAtomIndexSSC1,
                                         final int rootAtomIndexSSC2, final double shiftTol,
                                         final double matchFactorThrs) {
        //        System.out.println("\n"
        //                                   + ssc1.getHoseCodes()
        //                                         .get(rootAtomIndexSSC1));
        //        System.out.println(HOSECodeBuilder.buildConnectionTree(ssc1.getStructure(), rootAtomIndexSSC1, maxSphere));
        //        System.out.println("vs.");
        //        System.out.println(ssc2.getHoseCodes()
        //                               .get(rootAtomIndexSSC2));
        //        System.out.println(HOSECodeBuilder.buildConnectionTree(ssc2.getStructure(), rootAtomIndexSSC2, maxSphere));

        final List<SSC> extendedSSCList = new ArrayList<>();
        final int matchingSphere = AssemblyUtils.getMaximumMatchingSphereHOSECode(ssc1, ssc2, rootAtomIndexSSC1,
                                                                                  rootAtomIndexSSC2);
        final List<Integer> substructureAtomIndicesListHOSECodeSSC1 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc1.getStructure(), rootAtomIndexSSC1, matchingSphere);
        final List<Integer> substructureAtomIndicesListHOSECodeSSC2 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc2.getStructure(), rootAtomIndexSSC2, matchingSphere);
        final int matchingAtomCount = substructureAtomIndicesListHOSECodeSSC1.size();

        final Map<Integer, Integer> atomIndexMap = new HashMap<>();
        for (int i = 0; i
                < matchingAtomCount; i++) {
            atomIndexMap.put(substructureAtomIndicesListHOSECodeSSC1.get(i),
                             substructureAtomIndicesListHOSECodeSSC2.get(i));
        }
        //        System.out.println("\n -> overlapping sphere count: "
        //                                   + matchingSphere
        //                                   + " with "
        //                                   + matchingAtomCount
        //                                   + " atoms");
        //        System.out.println(" -> mapping: "
        //                                   + atomIndexMap);
        final ConnectionTree connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getStructure(),
                                                                                      rootAtomIndexSSC1, null);
        final ConnectionTree connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getStructure(),
                                                                                      rootAtomIndexSSC2, null);
        final boolean containsUnsaturatedAtomsSSC1 = connectionTreeSSC1.getNodeKeysInSphere(matchingSphere)
                                                                       .stream()
                                                                       .anyMatch(
                                                                               atomIndex -> ssc1.getUnsaturatedAtomIndices()
                                                                                                .contains(atomIndex));
        //        System.out.println("any unsaturated in sphere "
        //                                   + matchingSphere
        //                                   + " of SSC 1 to extend on? -> "
        //                                   + containsUnsaturatedAtomsSSC1);
        if (containsUnsaturatedAtomsSSC1) {
            int parentAtomIndexSSC2;
            try {
                final IAtomContainer extendedStructure = ssc1.getStructure()
                                                             .clone();
                final Spectrum extendedSpectrum = ssc1.getSpectrum()
                                                      .buildClone();
                final Assignment extendedAssignment = ssc1.getAssignment()
                                                          .buildClone();

                // add missing children and ring closures from connection tree in SSC2
                final List<Integer> addedAtomIndicesSSC1 = new ArrayList<>();
                for (final ConnectionTreeNode connectionTreeNodeSSC2 : connectionTreeSSC2.getNodesInSphere(
                        matchingSphere
                                + 1, true)) {
                    if (connectionTreeNodeSSC2.isRingClosureNode()) {
                        parentAtomIndexSSC2 = connectionTreeNodeSSC2.getParent()
                                                                    .getKey();
                        final int ringClosureParentAtomIndexSSC2 = connectionTreeNodeSSC2.getRingClosureParent()
                                                                                         .getKey();
                        final Integer parentAtomIndexSSC1 = utils.Utils.findKeyInMap(atomIndexMap, parentAtomIndexSSC2);
                        final Integer ringClosureParentAtomIndexSSC1 = utils.Utils.findKeyInMap(atomIndexMap,
                                                                                                ringClosureParentAtomIndexSSC2);
                        //                        System.out.println("--> ring closure to add: "
                        //                                                   + parentAtomIndexSSC1
                        //                                                   + " -> "
                        //                                                   + ringClosureParentAtomIndexSSC1);
                        if (parentAtomIndexSSC1
                                != null
                                && ringClosureParentAtomIndexSSC1
                                != null
                                && extendedStructure.getBond(extendedStructure.getAtom(parentAtomIndexSSC1),
                                                             extendedStructure.getAtom(ringClosureParentAtomIndexSSC1))
                                == null) {
                            //                            System.out.println(" --> adding ring closure allowed");
                            final IBond ringClosureBondSSC2 = ssc2.getStructure()
                                                                  .getBond(ssc2.getStructure()
                                                                               .getAtom(parentAtomIndexSSC2),
                                                                           ssc2.getStructure()
                                                                               .getAtom(ringClosureParentAtomIndexSSC2))
                                                                  .clone();
                            final IBond ringClosureBondSSC1 = new Bond(extendedStructure.getAtom(parentAtomIndexSSC1),
                                                                       extendedStructure.getAtom(
                                                                               ringClosureParentAtomIndexSSC1),
                                                                       ringClosureBondSSC2.getOrder());
                            ringClosureBondSSC1.setIsInRing(ringClosureBondSSC2.isInRing());
                            ringClosureBondSSC1.setIsAromatic(ringClosureBondSSC2.isAromatic());
                            if (Utils.isValidBondAddition(extendedStructure, parentAtomIndexSSC1,
                                                          ringClosureBondSSC1)) {
                                extendedStructure.addBond(ringClosureBondSSC1);
                            }
                        }
                    } else if (!atomIndexMap.containsValue(connectionTreeNodeSSC2.getKey())) {
                        final ConnectionTree subtreeToAddSSC2 = ConnectionTree.buildSubtree(connectionTreeSSC2,
                                                                                            connectionTreeNodeSSC2.getKey());
                        parentAtomIndexSSC2 = connectionTreeNodeSSC2.getParent()
                                                                    .getKey();
                        final Integer parentAtomIndexSSC1 = utils.Utils.findKeyInMap(atomIndexMap, parentAtomIndexSSC2);
                        final IBond bondToLink = ssc2.getStructure()
                                                     .getBond(ssc2.getStructure()
                                                                  .getAtom(parentAtomIndexSSC2), ssc2.getStructure()
                                                                                                     .getAtom(
                                                                                                             connectionTreeNodeSSC2.getKey()))
                                                     .clone();
                        //                        System.out.println("--> to add node: "
                        //                                                   + connectionTreeNodeSSC2.getKey()
                        //                                                   + "\n--->"
                        //                                                   + subtreeToAddSSC2);
                        //                        System.out.println("---> atom count before: "
                        //                                                   + extendedStructure.getAtomCount());
                        if (Utils.isValidBondAddition(extendedStructure, parentAtomIndexSSC1, bondToLink)) {
                            HOSECodeBuilder.addToAtomContainer(subtreeToAddSSC2, extendedStructure, parentAtomIndexSSC1,
                                                               bondToLink);
                            //                        System.out.println("---> atom count after: "
                            //                                                   + extendedStructure.getAtomCount());
                            // add indices to index map
                            int atomCounterSSC1 = extendedStructure.getAtomCount()
                                    - subtreeToAddSSC2.getNodesCount(false);
                            for (int s = 0; s
                                    <= subtreeToAddSSC2.getMaxSphere(); s++) {
                                // first add all atoms and its parents (previous sphere only, incl. bonds) to structure
                                final List<ConnectionTreeNode> nodesInSphere = subtreeToAddSSC2.getNodesInSphere(s,
                                                                                                                 false);
                                for (int k = 0; k
                                        < nodesInSphere.size(); k++) {
                                    atomIndexMap.put(atomCounterSSC1, nodesInSphere.get(k)
                                                                                   .getKey());
                                    addedAtomIndicesSSC1.add(atomCounterSSC1);
                                    final Signal signalToAddSSC2 = ssc2.getSpectrum()
                                                                       .getSignal(ssc2.getAssignment()
                                                                                      .getIndex(0, nodesInSphere.get(k)
                                                                                                                .getKey()));
                                    if (signalToAddSSC2
                                            != null) {
                                        AssemblyUtils.addSignalToSSC(extendedSpectrum, extendedAssignment,
                                                                     signalToAddSSC2.buildClone(), atomCounterSSC1);
                                    }

                                    atomCounterSSC1++;
                                }
                            }
                        }
                    }
                }
                //                System.out.println(" -> mapping 2: "
                //                                           + atomIndexMap);

                // add possible missing bonds between added atoms from tree and already existing atoms in structure
                IAtom mappedAtomSSC2;
                IBond bond, bondToAdd;
                for (final int mappedAtomIndexSSC2 : atomIndexMap.values()) {
                    mappedAtomSSC2 = ssc2.getStructure()
                                         .getAtom(mappedAtomIndexSSC2);
                    for (final IAtom mappedAtomNeighborSSC2 : ssc2.getStructure()
                                                                  .getConnectedAtomsList(mappedAtomSSC2)) {
                        if (atomIndexMap.containsValue(mappedAtomNeighborSSC2.getIndex())) {
                            final int key1 = utils.Utils.findKeyInMap(atomIndexMap, mappedAtomIndexSSC2);
                            final int key2 = utils.Utils.findKeyInMap(atomIndexMap, mappedAtomNeighborSSC2.getIndex());
                            if (extendedStructure.getBond(extendedStructure.getAtom(key1),
                                                          extendedStructure.getAtom(key2))
                                    == null) {
                                bond = ssc2.getStructure()
                                           .getBond(mappedAtomSSC2, mappedAtomNeighborSSC2)
                                           .clone();
                                bondToAdd = new Bond(extendedStructure.getAtom(key1), extendedStructure.getAtom(key2),
                                                     bond.getOrder());
                                bondToAdd.setIsInRing(bond.isInRing());
                                bondToAdd.setIsAromatic(bond.isAromatic());
                                extendedStructure.addBond(bondToAdd);
                            }
                        }
                    }
                }
                final SSC extendedSSC = new SSC();
                extendedSSC.setStructure(extendedStructure);
                extendedSSC.setSpectrum(extendedSpectrum);
                extendedSSC.setAssignment(extendedAssignment);
                extendedSSC.setUnsaturatedAtomIndices(
                        utils.Utils.getUnsaturatedAtomIndices(extendedSSC.getStructure()));
                extendedSSC.setHoseCodes(utils.Utils.buildHOSECodes(extendedSSC.getStructure(), maxSphere));
                final boolean isValidExtension = AssemblyUtils.isValidExtension(extendedSSC, querySpectrum, shiftTol,
                                                                                matchFactorThrs, mf,
                                                                                ssc1.getStructure());
                //                System.out.println("isValidExtension ? -> "
                //                                           + isValidExtension);
                if (isValidExtension) {
                    extendedSSCList.add(extendedSSC);
                }
            } catch (final CloneNotSupportedException e) {
                e.printStackTrace();
            }
        }

        return extendedSSCList;
    }

    public static Map<String, SSC> assembleWithBacktracking(final List<SSC> sscList, final Spectrum querySpectrum,
                                                            final String mf, final double shiftTol,
                                                            final double matchFactorThrs, final int minMatchingSphere,
                                                            final int nStarts) {
        final Map<String, SSC> solutions = new HashMap<>();
        Stack<Object[]> intermediates;
        Object[] stackItem;
        SSC intermediate;
        int currentSSCIndex;
        List<SSC> extendedSSCList;
        String smiles;
        Set<String> smilesSet;
        boolean stop;
        for (int i = 0; i
                < sscList.size(); i++) {
            if (i
                    >= nStarts) {
                break;
            }

            if (AssemblyUtils.isFinalSSC(sscList.get(i), querySpectrum, shiftTol, matchFactorThrs, mf)) {
                try {
                    smiles = smilesGenerator.create(sscList.get(i)
                                                           .getStructure());
                    if (!solutions.containsKey(smiles)) {
                        try {
                            depictionGenerator.depict(sscList.get(i)
                                                             .getStructure())
                                              .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                               + i
                                                               + "-"
                                                               + "X"
                                                               + "_"
                                                               + smiles
                                                               + ".png");
                        } catch (final IOException | CDKException e) {
                            e.printStackTrace();
                        }
                    }

                    solutions.putIfAbsent(smiles, sscList.get(i));
                    break;
                } catch (final CDKException e) {
                    e.printStackTrace();
                }
            }

            intermediates = new Stack<>();
            intermediates.push(new Object[]{sscList.get(i).buildClone(), 0, new HashSet<String>()});

            while (!intermediates.empty()) {
                stackItem = intermediates.pop();
                intermediate = (SSC) stackItem[0];
                currentSSCIndex = (int) stackItem[1];
                smilesSet = (Set<String>) stackItem[2];

                if (currentSSCIndex
                        >= sscList.size()) {
                    continue;
                }
                stop = false;
                for (int j = currentSSCIndex; j
                        < sscList.size(); j++) {
                    System.out.println("\n\n--> ssc pair: "
                                               + i
                                               + " vs. "
                                               + j);
                    extendedSSCList = Assembly.assemblyOverlaps(intermediate, sscList.get(j), querySpectrum, mf, null,
                                                                shiftTol, matchFactorThrs, minMatchingSphere);

                    AssemblyUtils.sortExtendedSSCList(extendedSSCList);
                    // remove duplicates after sorting since we do not check the spectral similarities when removing them
                    utils.Utils.removeDuplicatesFromSSCList(extendedSSCList, smilesSet);//new HashSet<>());
                    //                    Collections.reverse(extendedSSCList);
                    if (!extendedSSCList.isEmpty()) {
                        //                        for (final SSC extendedSSC : extendedSSCList) {
                        final SSC extendedSSC = extendedSSCList.get(0);
                        if (AssemblyUtils.isFinalSSC(extendedSSC, querySpectrum, shiftTol, matchFactorThrs, mf)) {
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                if (!solutions.containsKey(smiles)) {
                                    try {
                                        depictionGenerator.depict(extendedSSC.getStructure())
                                                          .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                                           + i
                                                                           + "-"
                                                                           + j
                                                                           + "_"
                                                                           + smiles
                                                                           + ".png");
                                    } catch (final IOException | CDKException e) {
                                        e.printStackTrace();
                                    }
                                }
                                solutions.putIfAbsent(smiles, extendedSSC);
                                stop = true;
                            } catch (final CDKException e) {
                                e.printStackTrace();
                            }
                        } else {
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                depictionGenerator.depict(extendedSSC.getStructure())
                                                  .writeTo("/Users/mwenk/Downloads/depictions/extended_"
                                                                   + i
                                                                   + "-"
                                                                   + j
                                                                   + "_"
                                                                   + smiles
                                                                   + "_"
                                                                   + extendedSSC.getUnsaturatedAtomIndices()
                                                                   + ".png");
                            } catch (final IOException | CDKException e) {
                                e.printStackTrace();
                            }

                            intermediates.push(new Object[]{extendedSSC.buildClone(), j
                                    + 1, new HashSet<>(smilesSet)});
                        }
                        //                        }
                    }
                    if (stop) {
                        break;
                    }
                }
                if (stop) {
                    break;
                }
            }
        }

        return solutions;
    }

    public static Map<String, SSC> assembleSequentially(final List<SSC> sscList, final Spectrum querySpectrum,
                                                        final String mf, final double shiftTol,
                                                        final double matchFactorThrs, final int minMatchingSphere,
                                                        final int nStarts) {
        final Map<String, SSC> solutions = new HashMap<>();

        List<SSC> extendedSSCList;
        SSC intermediate, extendedSSC;
        String smiles;
        boolean stop = false;
        for (int i = 0; i
                < sscList.size(); i++) {
            if (i
                    >= nStarts) {
                break;
            }
            intermediate = sscList.get(i)
                                  .buildClone();
            if (AssemblyUtils.isFinalSSC(intermediate, querySpectrum, shiftTol, matchFactorThrs, mf)) {
                try {
                    smiles = smilesGenerator.create(intermediate.getStructure());
                    if (!solutions.containsKey(smiles)) {
                        try {
                            depictionGenerator.depict(intermediate.getStructure())
                                              .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                               + i
                                                               + "-"
                                                               + "X"
                                                               + "_"
                                                               + smiles
                                                               + ".png");
                        } catch (final IOException | CDKException e) {
                            e.printStackTrace();
                        }
                    }

                    solutions.putIfAbsent(smiles, intermediate);
                    stop = true;
                } catch (final CDKException e) {
                    e.printStackTrace();
                }
            } else {
                for (int j = //0
                     //
                     i
                             + 1
                     //
                     ; j
                             < sscList.size(); j++) {
                    if (j
                            == i) {
                        continue;
                    }
                    extendedSSCList = Assembly.assemblyOverlaps(intermediate, sscList.get(j), querySpectrum, mf, null,
                                                                shiftTol, matchFactorThrs, minMatchingSphere);
                    AssemblyUtils.sortExtendedSSCList(extendedSSCList);
                    if (!extendedSSCList.isEmpty()) {
                        extendedSSC = extendedSSCList.get(0);
                        if (AssemblyUtils.isFinalSSC(extendedSSC, querySpectrum, shiftTol, matchFactorThrs, mf)) {
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                if (!solutions.containsKey(smiles)) {
                                    try {
                                        depictionGenerator.depict(extendedSSC.getStructure())
                                                          .writeTo("/Users/mwenk/Downloads/depictions/final_"
                                                                           + i
                                                                           + "-"
                                                                           + j
                                                                           + "_"
                                                                           + smiles
                                                                           + ".png");
                                    } catch (final IOException | CDKException e) {
                                        e.printStackTrace();
                                    }
                                }

                                solutions.putIfAbsent(smiles, extendedSSC);
                                stop = true;
                            } catch (final CDKException e) {
                                e.printStackTrace();
                            }
                        } else {
                            try {
                                smiles = smilesGenerator.create(extendedSSC.getStructure());
                                depictionGenerator.depict(extendedSSC.getStructure())
                                                  .writeTo("/Users/mwenk/Downloads/depictions/extended_"
                                                                   + i
                                                                   + "-"
                                                                   + j
                                                                   + "_"
                                                                   + smiles
                                                                   + "_"
                                                                   + extendedSSC.getUnsaturatedAtomIndices()
                                                                   + ".png");
                            } catch (final IOException | CDKException e) {
                                e.printStackTrace();
                            }
                            if (!intermediate.getUnsaturatedAtomIndices()
                                             .isEmpty()) {
                                intermediate = extendedSSC.buildClone();
                            }
                        }
                    }
                    if (stop) {
                        break;
                    }
                }
            }
            if (stop) {
                break;
            }
        }

        return solutions;
    }
}
