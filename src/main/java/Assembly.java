import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.hose.Utils;
import casekit.nmr.hose.model.ConnectionTree;
import casekit.nmr.hose.model.ConnectionTreeNode;
import casekit.nmr.model.Assignment;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import model.SSC;
import org.openscience.cdk.Bond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import utils.AssemblyUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Assembly {

    public static List<SSC> assemblyOverlaps(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum,
                                             final String mf, final Integer maxSphere, final double shiftTol,
                                             final double matchFactorThrs, final int minMatchingSphereCount) {
        final Map<Integer, List<Integer[]>> overlaps = AssemblyUtils.getOverlaps(ssc1, ssc2, minMatchingSphereCount);
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
        final List<Integer> substructureAtomIndicesListHOSECodeSSC1 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc1.getStructure(), rootAtomIndexSSC1, maxSphere);
        //        System.out.println(substructureAtomIndicesListHOSECodeSSC1);
        final List<String> hoseCodeSpheresSSC1 = Utils.splitHOSECodeIntoSpheres(ssc1.getHoseCodes()
                                                                                    .get(rootAtomIndexSSC1));
        final List<Integer> substructureAtomIndicesListHOSECodeSSC2 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc2.getStructure(), rootAtomIndexSSC2, maxSphere);
        //        System.out.println(substructureAtomIndicesListHOSECodeSSC2);
        final List<String> hoseCodeSpheresSSC2 = Utils.splitHOSECodeIntoSpheres(ssc2.getHoseCodes()
                                                                                    .get(rootAtomIndexSSC2));
        int matchingSphereCount = 0;
        int matchingAtomCount = 0;
        for (int i = 0; i
                < Math.min(hoseCodeSpheresSSC1.size(), hoseCodeSpheresSSC2.size()); i++) {
            if (!hoseCodeSpheresSSC1.get(i)
                                    .equals(hoseCodeSpheresSSC2.get(i))) {
                break;
            }
            matchingSphereCount++;
        }
        for (int s = 0; s
                < matchingSphereCount; s++) {
            matchingAtomCount += Utils.countAtoms(hoseCodeSpheresSSC1.get(s));
        }

        final Map<Integer, Integer> atomIndexMap = new HashMap<>();
        for (int i = 0; i
                < matchingAtomCount; i++) {
            atomIndexMap.put(substructureAtomIndicesListHOSECodeSSC1.get(i),
                             substructureAtomIndicesListHOSECodeSSC2.get(i));
        }
        //        System.out.println("\n -> overlapping sphere count: "
        //                                   + matchingSphereCount
        //                                   + " with "
        //                                   + matchingAtomCount
        //                                   + " atoms");
        //        System.out.println(" -> mapping: "
        //                                   + atomIndexMap);
        final boolean containsUnsaturatedAtomsSSC1 = atomIndexMap.keySet()
                                                                 .stream()
                                                                 .anyMatch(atomIndex -> ssc1.getUnsaturatedAtomIndices()
                                                                                            .contains(atomIndex));
        //        System.out.println("any unsaturated in SSC to extend? -> "
        //                                   + containsUnsaturatedAtomsSSC1);
        //        System.out.println(matchingSphereCount
        //                                   + " < "
        //                                   + hoseCodeSpheresSSC1.size()
        //                                   + " -> "
        //                                   + (matchingSphereCount
        //                < hoseCodeSpheresSSC1.size()));

        if (containsUnsaturatedAtomsSSC1
                && matchingSphereCount
                < hoseCodeSpheresSSC1.size()
                && matchingSphereCount
                < hoseCodeSpheresSSC2.size()) {
            // extend atom mapping within next sphere (from maximum matching sphere)
            final Map<Integer, ArrayList<String>> splitPositionsSSC1 = Utils.splitHOSECodeSphereIntoPositions(
                    hoseCodeSpheresSSC1.get(matchingSphereCount), false);
            final Map<Integer, ArrayList<String>> splitPositionsSSC2 = Utils.splitHOSECodeSphereIntoPositions(
                    hoseCodeSpheresSSC2.get(matchingSphereCount), false);
            int atomIndexInHOSECodeSSC1 = matchingAtomCount;
            int atomIndexInHOSECodeSSC2 = matchingAtomCount;
            int elemCountPrevSphere, parentAtomIndexSSC2;
            boolean isEmptySSC1, stop;
            Map<Integer, ArrayList<String>> splitPositionsPrevSpheresSSC2;
            List<String> positionElementsSSC2;
            for (final int posIndexSSC2 : splitPositionsSSC2.keySet()) {
                // ######## detection of parent atom in SSC2
                for (int k = 0; k
                        < matchingSphereCount; k++) {
                    splitPositionsPrevSpheresSSC2 = Utils.splitHOSECodeSphereIntoPositions(hoseCodeSpheresSSC2.get(k), k
                            == 0);
                    if (k
                            >= matchingSphereCount
                            - 1) {
                        elemCountPrevSphere = 0;
                        stop = false;
                        for (final int posIndexPrevSphere : splitPositionsPrevSpheresSSC2.keySet()) {
                            positionElementsSSC2 = splitPositionsPrevSpheresSSC2.get(posIndexPrevSphere);
                            for (int i = 0; i
                                    < positionElementsSSC2.size(); i++) {
                                if (elemCountPrevSphere
                                        == posIndexSSC2) {
                                    stop = true;
                                    break;
                                }
                                elemCountPrevSphere++;
                            }
                            if (stop) {
                                break;
                            }
                        }
                    }
                }
                // ########
                positionElementsSSC2 = splitPositionsSSC2.get(posIndexSSC2);
                for (int elemIndexSSC2 = 0; elemIndexSSC2
                        < positionElementsSSC2.size(); elemIndexSSC2++) {
                    isEmptySSC1 = false;
                    if (!splitPositionsSSC1.containsKey(posIndexSSC2)
                            || elemIndexSSC2
                            >= splitPositionsSSC1.get(posIndexSSC2)
                                                 .size()
                            || splitPositionsSSC1.get(posIndexSSC2)
                                                 .get(elemIndexSSC2)
                            == null
                            || elemIndexSSC2
                            >= splitPositionsSSC1.get(posIndexSSC2)
                                                 .size()
                            || Utils.countAtoms(splitPositionsSSC1.get(posIndexSSC2)
                                                                  .get(elemIndexSSC2))
                            == 0) {
                        isEmptySSC1 = true;
                    }
                    if (!isEmptySSC1) {
                        if (atomIndexInHOSECodeSSC1
                                < substructureAtomIndicesListHOSECodeSSC1.size()
                                && atomIndexInHOSECodeSSC2
                                < substructureAtomIndicesListHOSECodeSSC2.size()) {
                            atomIndexMap.put(substructureAtomIndicesListHOSECodeSSC1.get(atomIndexInHOSECodeSSC1),
                                             substructureAtomIndicesListHOSECodeSSC2.get(atomIndexInHOSECodeSSC2));
                        }
                        atomIndexInHOSECodeSSC1++;
                        atomIndexInHOSECodeSSC2++;
                    } else if (positionElementsSSC2.get(elemIndexSSC2)
                            != null
                            && Utils.countAtoms(positionElementsSSC2.get(elemIndexSSC2))
                            > 0) {
                        atomIndexInHOSECodeSSC2++;
                    }
                }
            }
            //            System.out.println(" -> mapping 2: "
            //                                       + atomIndexMap);

            try {
                final IAtomContainer extendedStructure = ssc1.getStructure()
                                                             .clone();
                final Spectrum extendedSpectrum = ssc1.getSpectrum()
                                                      .buildClone();
                final Assignment extendedAssignment = ssc1.getAssignment()
                                                          .buildClone();

                // add missing children and ring closures from connection tree in SSC2
                final ConnectionTree connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getStructure(),
                                                                                              rootAtomIndexSSC2,
                                                                                              maxSphere);
                for (final ConnectionTreeNode connectionTreeNodeSSC2 : connectionTreeSSC2.getNodes(true)) {
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
                            extendedStructure.addBond(ringClosureBondSSC1);
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
                            final List<ConnectionTreeNode> nodesInSphere = subtreeToAddSSC2.getNodesInSphere(s, false);
                            for (int k = 0; k
                                    < nodesInSphere.size(); k++) {
                                atomIndexMap.put(atomCounterSSC1, nodesInSphere.get(k)
                                                                               .getKey());
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
                //                System.out.println(" -> mapping 3: "
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
}
