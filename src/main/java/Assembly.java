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

import java.util.*;

public class Assembly {

    public static List<SSC> assemblyOverlaps(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum,
                                             final String mf, final int maxSphere, final double shiftTol,
                                             final double matchFactorThrs, final int minMatchingSphereCount) {
        final Map<Integer, List<Integer[]>> overlaps = AssemblyUtils.getOverlaps(ssc1, ssc2, minMatchingSphereCount);
        final List<SSC> extendedSSCList = new ArrayList<>();
        for (final Map.Entry<Integer, List<Integer[]>> integerListEntry : overlaps.entrySet()) {
            System.out.println("---> sphere match: "
                                       + integerListEntry.getKey());
            for (final Integer[] rootAtomIndices : integerListEntry.getValue()) {
                System.out.println("\n--> root atom pair: "
                                           + Arrays.toString(rootAtomIndices));
                extendedSSCList.addAll(
                        Assembly.assemblyCore(querySpectrum, mf, maxSphere, ssc1, ssc2, rootAtomIndices[0],
                                              rootAtomIndices[1], shiftTol, matchFactorThrs));
            }
        }

        return extendedSSCList;
    }

    public static List<SSC> assemblyCore(final Spectrum querySpectrum, final String mf, final int maxSphere,
                                         final SSC ssc1, final SSC ssc2, final int rootAtomIndexSSC1,
                                         final int rootAtomIndexSSC2, final double shiftTol,
                                         final double matchFactorThrs) {
        System.out.println(ssc1.getHoseCodes()
                               .get(rootAtomIndexSSC1));
        System.out.println("vs.");
        System.out.println(ssc2.getHoseCodes()
                               .get(rootAtomIndexSSC2));
        final List<SSC> extendedSSCList = new ArrayList<>();
        final List<Integer> substructureAtomIndicesListHOSECodeSSC1 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc1.getStructure(), rootAtomIndexSSC1, maxSphere);
        final List<String> hoseCodeSpheresSSC1 = Utils.splitHOSECodeIntoSpheres(ssc1.getHoseCodes()
                                                                                    .get(rootAtomIndexSSC1));
        final List<Integer> substructureAtomIndicesListHOSECodeSSC2 = Fragmentation.buildSubstructureAtomIndicesList(
                ssc2.getStructure(), rootAtomIndexSSC2, maxSphere);
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
        System.out.println("\n -> overlapping sphere count: "
                                   + matchingSphereCount
                                   + " with "
                                   + matchingAtomCount
                                   + " atoms");
        System.out.println(" -> mapping: "
                                   + atomIndexMap);
        final boolean containsUnsaturatedAtomsSSC1 = atomIndexMap.keySet()
                                                                 .stream()
                                                                 .anyMatch(atomIndex -> ssc1.getUnsaturatedAtomIndices()
                                                                                            .contains(atomIndex));
        System.out.println("any unsaturated in SSC to extend? -> "
                                   + containsUnsaturatedAtomsSSC1);
        System.out.println(matchingSphereCount
                                   + " < "
                                   + hoseCodeSpheresSSC1.size()
                                   + " -> "
                                   + (matchingSphereCount
                < hoseCodeSpheresSSC1.size()));

        if (containsUnsaturatedAtomsSSC1
                && matchingSphereCount
                < hoseCodeSpheresSSC1.size()
                && matchingSphereCount
                < hoseCodeSpheresSSC2.size()) {
            final Map<Integer, ArrayList<String>> splitPositionsSSC1 = Utils.splitHOSECodeSphereIntoPositions(
                    hoseCodeSpheresSSC1.get(matchingSphereCount), false);
            final Map<Integer, ArrayList<String>> splitPositionsSSC2 = Utils.splitHOSECodeSphereIntoPositions(
                    hoseCodeSpheresSSC2.get(matchingSphereCount), false);
            int atomIndexSSC1 = matchingAtomCount;
            int atomIndexSSC2 = matchingAtomCount;
            int totalElemCount, elemCountPrevSphere, parentAtomIndexSSC2;
            boolean isEmptySSC1, stop;
            Map<Integer, ArrayList<String>> splitPositionsPrevSpheresSSC2;
            List<String> positionElementsSSC2;
            final Map<Integer, List<Integer>> atomIndicesToAddSSC2 = new HashMap<>();
            for (final int posIndexSSC2 : splitPositionsSSC2.keySet()) {
                // ######## detection of parent atom in SSC2
                totalElemCount = 0;
                parentAtomIndexSSC2 = -1;
                for (int k = 0; k
                        < matchingSphereCount; k++) {
                    splitPositionsPrevSpheresSSC2 = Utils.splitHOSECodeSphereIntoPositions(hoseCodeSpheresSSC2.get(k), k
                            == 0);
                    if (k
                            < matchingSphereCount
                            - 1) {
                        for (final int posIndexPrevSphere : splitPositionsPrevSpheresSSC2.keySet()) {
                            totalElemCount += Utils.countAtoms(splitPositionsPrevSpheresSSC2.get(posIndexPrevSphere)
                                                                                            .stream()
                                                                                            .reduce("",
                                                                                                    (total, pos) -> total
                                                                                                            + pos));
                        }
                    } else {
                        elemCountPrevSphere = 0;
                        stop = false;
                        for (final int posIndexPrevSphere : splitPositionsPrevSpheresSSC2.keySet()) {
                            positionElementsSSC2 = splitPositionsPrevSpheresSSC2.get(posIndexPrevSphere);
                            for (int i = 0; i
                                    < positionElementsSSC2.size(); i++) {
                                if (elemCountPrevSphere
                                        == posIndexSSC2) {
                                    parentAtomIndexSSC2 = substructureAtomIndicesListHOSECodeSSC2.get(totalElemCount
                                                                                                              + elemCountPrevSphere);
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
                        if (atomIndexSSC1
                                < substructureAtomIndicesListHOSECodeSSC1.size()
                                && atomIndexSSC2
                                < substructureAtomIndicesListHOSECodeSSC2.size()) {
                            atomIndexMap.put(substructureAtomIndicesListHOSECodeSSC1.get(atomIndexSSC1),
                                             substructureAtomIndicesListHOSECodeSSC2.get(atomIndexSSC2));
                        }
                        atomIndexSSC1++;
                        atomIndexSSC2++;
                    } else if (positionElementsSSC2.get(elemIndexSSC2)
                            != null
                            && Utils.countAtoms(positionElementsSSC2.get(elemIndexSSC2))
                            > 0) {
                        if (parentAtomIndexSSC2
                                >= 0
                                && atomIndexSSC2
                                < substructureAtomIndicesListHOSECodeSSC2.size()) {
                            atomIndicesToAddSSC2.putIfAbsent(parentAtomIndexSSC2, new ArrayList<>());
                            atomIndicesToAddSSC2.get(parentAtomIndexSSC2)
                                                .add(substructureAtomIndicesListHOSECodeSSC2.get(atomIndexSSC2));
                        }
                        atomIndexSSC2++;
                    }
                }
            }
            System.out.println(" -> mapping: "
                                       + atomIndexMap);
            boolean allBondAdditionsAreValid = true;
            try {
                final IAtomContainer extendedStructure = ssc1.getStructure()
                                                             .clone();
                final Spectrum extendedSpectrum = ssc1.getSpectrum()
                                                      .buildClone();
                final Assignment extendedAssignment = ssc1.getAssignment()
                                                          .buildClone();
                for (final int parentAtomIndexSSC2Temp : atomIndicesToAddSSC2.keySet()) {
                    for (int i = 0; i
                            < atomIndicesToAddSSC2.get(parentAtomIndexSSC2Temp)
                                                  .size(); i++) {
                        final int parentAtomIndexSSC1 = utils.Utils.findKeyInMap(atomIndexMap, parentAtomIndexSSC2Temp);
                        final int atomIndexToAddSSC2 = atomIndicesToAddSSC2.get(parentAtomIndexSSC2Temp)
                                                                           .get(i);
                        IBond bond = ssc2.getStructure()
                                         .getBond(ssc2.getStructure()
                                                      .getAtom(parentAtomIndexSSC2Temp), ssc2.getStructure()
                                                                                             .getAtom(
                                                                                                     atomIndexToAddSSC2));
                        boolean isValidBondAddition = bond
                                != null
                                && casekit.nmr.Utils.isValidBondAddition(ssc1.getStructure(), parentAtomIndexSSC1,
                                                                         bond);
                        if (!isValidBondAddition) {
                            allBondAdditionsAreValid = false;
                            break;
                        }
                        final IAtom atomToAdd = ssc2.getStructure()
                                                    .getAtom(atomIndexToAddSSC2)
                                                    .clone();
                        IBond bondToAdd = new Bond(extendedStructure.getAtom(parentAtomIndexSSC1), atomToAdd,
                                                   bond.getOrder());
                        bondToAdd.setIsInRing(bond.isInRing());
                        bondToAdd.setIsAromatic(bond.isAromatic());
                        extendedStructure.addAtom(atomToAdd);
                        extendedStructure.addBond(bondToAdd);
                        final int indexOfAddedAtomSSC1 = extendedStructure.getAtomCount()
                                - 1;
                        Signal signalToAddSSC2 = ssc2.getSpectrum()
                                                     .getSignal(ssc2.getAssignment()
                                                                    .getIndex(0, atomIndexToAddSSC2));
                        if (signalToAddSSC2
                                != null) {
                            AssemblyUtils.addSignalToSSC(extendedSpectrum, extendedAssignment, signalToAddSSC2,
                                                         indexOfAddedAtomSSC1);
                        }
                        atomIndexMap.put(indexOfAddedAtomSSC1, atomIndexToAddSSC2);
                        System.out.println("atomIndexMap new: "
                                                   + atomIndexMap);

                        // check for ring closure at added atom from SSC2
                        final IAtom addedAtomInSSC2 = ssc2.getStructure()
                                                          .getAtom(atomIndexToAddSSC2);
                        if (addedAtomInSSC2.isInRing()) {
                            for (final IBond bondInSSC2 : addedAtomInSSC2.bonds()) {
                                final Integer keyBondAtom1 = utils.Utils.findKeyInMap(atomIndexMap,
                                                                                      bondInSSC2.getAtom(0)
                                                                                                .getIndex());
                                final Integer keyBondAtom2 = utils.Utils.findKeyInMap(atomIndexMap,
                                                                                      bondInSSC2.getAtom(1)
                                                                                                .getIndex());
                                if (keyBondAtom1
                                        != null
                                        && keyBondAtom2
                                        != null) {
                                    // check whether bond does not exist
                                    if (extendedStructure.getBond(extendedStructure.getAtom(keyBondAtom1),
                                                                  extendedStructure.getAtom(keyBondAtom2))
                                            == null) {
                                        final IBond ringClosureBond = new Bond(extendedStructure.getAtom(keyBondAtom1),
                                                                               extendedStructure.getAtom(keyBondAtom2),
                                                                               bondInSSC2.getOrder());
                                        ringClosureBond.setIsInRing(bondInSSC2.isInRing());
                                        ringClosureBond.setIsAromatic(bondInSSC2.isAromatic());

                                        if (casekit.nmr.Utils.isValidBondAddition(extendedStructure, keyBondAtom1,
                                                                                  ringClosureBond)
                                                && casekit.nmr.Utils.isValidBondAddition(extendedStructure,
                                                                                         keyBondAtom2,
                                                                                         ringClosureBond)) {
                                            extendedStructure.addBond(ringClosureBond);
                                        }
                                    }
                                }
                            }
                        }

                        // for each neighbor of added atom build a connection tree with still unvisited nodes to add
                        for (final IAtom connectedAtomSSC2 : ssc2.getStructure()
                                                                 .getConnectedAtomsList(addedAtomInSSC2)) {
                            if (!atomIndexMap.containsValue(connectedAtomSSC2.getIndex())) {
                                // try to append the rest of connected but not visited atoms in SSC2
                                final ConnectionTree connectionTreeToAddSSC2 = HOSECodeBuilder.buildConnectionTree(
                                        ssc2.getStructure(), connectedAtomSSC2.getIndex(), null,
                                        new HashSet<>(atomIndexMap.values()));
                                // add atoms from tree to container
                                HOSECodeBuilder.addToAtomContainer(connectionTreeToAddSSC2, extendedStructure,
                                                                   indexOfAddedAtomSSC1, ssc2.getStructure()
                                                                                             .getBond(addedAtomInSSC2,
                                                                                                      connectedAtomSSC2)
                                                                                             .clone());
                                // add indices to index map
                                int atomCounterSSC1 = extendedStructure.getAtomCount()
                                        - connectionTreeToAddSSC2.getNodesCount(false);
                                for (int s = 0; s
                                        <= connectionTreeToAddSSC2.getMaxSphere(); s++) {
                                    // first add all atoms and its parents (previous sphere only, incl. bonds) to structure
                                    final List<ConnectionTreeNode> nodesInSphere = connectionTreeToAddSSC2.getNodesInSphere(
                                            s, false);
                                    for (int k = 0; k
                                            < nodesInSphere.size(); k++) {
                                        atomIndexMap.put(atomCounterSSC1, nodesInSphere.get(k)
                                                                                       .getKey());
                                        signalToAddSSC2 = ssc2.getSpectrum()
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
                        System.out.println("atomIndexMap final: "
                                                   + atomIndexMap);
                        // add possible missing bonds between added atoms from tree and already existing atoms in structure
                        IAtom mappedAtomSSC2;
                        for (final int mappedAtomIndexSSC2 : atomIndexMap.values()) {
                            mappedAtomSSC2 = ssc2.getStructure()
                                                 .getAtom(mappedAtomIndexSSC2);
                            for (final IAtom mappedAtomNeighborSSC2 : ssc2.getStructure()
                                                                          .getConnectedAtomsList(mappedAtomSSC2)) {
                                if (atomIndexMap.containsValue(mappedAtomNeighborSSC2.getIndex())) {
                                    final int key1 = utils.Utils.findKeyInMap(atomIndexMap, mappedAtomIndexSSC2);
                                    final int key2 = utils.Utils.findKeyInMap(atomIndexMap,
                                                                              mappedAtomNeighborSSC2.getIndex());
                                    if (extendedStructure.getBond(extendedStructure.getAtom(key1),
                                                                  extendedStructure.getAtom(key2))
                                            == null) {
                                        bond = ssc2.getStructure()
                                                   .getBond(mappedAtomSSC2, mappedAtomNeighborSSC2)
                                                   .clone();
                                        isValidBondAddition = casekit.nmr.Utils.isValidBondAddition(extendedStructure,
                                                                                                    parentAtomIndexSSC1,
                                                                                                    bond);
                                        if (isValidBondAddition) {
                                            bondToAdd = new Bond(extendedStructure.getAtom(key1),
                                                                 extendedStructure.getAtom(key2), bond.getOrder());
                                            bondToAdd.setIsInRing(bond.isInRing());
                                            bondToAdd.setIsAromatic(bond.isAromatic());

                                            extendedStructure.addBond(bondToAdd);
                                        } else {
                                            allBondAdditionsAreValid = false;
                                            break;
                                        }
                                    }
                                }
                            }
                            if (!allBondAdditionsAreValid) {
                                break;
                            }
                        }
                    }
                    if (!allBondAdditionsAreValid) {
                        break;
                    }
                }
                System.out.println("areAllBondAdditionValid ? -> "
                                           + allBondAdditionsAreValid);
                if (allBondAdditionsAreValid) {
                    System.out.println("extended spectrum -> "
                                               + extendedSpectrum);

                    final SSC extendedSSC = new SSC();
                    extendedSSC.setStructure(extendedStructure);
                    extendedSSC.setSpectrum(extendedSpectrum);
                    extendedSSC.setAssignment(extendedAssignment);
                    extendedSSC.setUnsaturatedAtomIndices(
                            utils.Utils.getUnsaturatedAtomIndices(extendedSSC.getStructure()));
                    extendedSSC.setHoseCodes(utils.Utils.buildHOSECodes(extendedSSC.getStructure(), maxSphere));
                    final boolean isValidExtension = AssemblyUtils.isValidExtension(extendedSSC, querySpectrum,
                                                                                    shiftTol, matchFactorThrs, mf,
                                                                                    ssc1.getStructure());
                    System.out.println("isValidExtension ? -> "
                                               + isValidExtension);
                    if (isValidExtension) {
                        extendedSSCList.add(extendedSSC);
                    }
                }
            } catch (final CloneNotSupportedException e) {
                e.printStackTrace();
            }

        }

        return extendedSSCList;
    }
}
