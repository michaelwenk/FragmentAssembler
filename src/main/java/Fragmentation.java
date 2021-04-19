import casekit.nmr.model.Assignment;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Utils;
import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import utils.AssemblyUtils;
import utils.ParallelTasks;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

public class Fragmentation {

    /**
     * Builds a collection of substructure-subspectrum-correlations (SSC objects) from
     * an atom container set for all its molecules and atoms by using a
     * breadth first search with spherical limit.
     *
     * @param dataSets  Set of components for each structure to build SSCs from it
     * @param maxSphere Spherical limit for building a substructure into
     *                  all directions
     * @param nThreads  Number of threads to use for parallelization
     *
     * @return
     *
     * @throws java.lang.InterruptedException
     * @see Fragmentation#buildSSCCollection(IAtomContainer, Spectrum, Assignment, int)
     */
    public static List<SSC> buildSSCCollection(final List<DataSet> dataSets, final int maxSphere,
                                               final int nThreads) throws InterruptedException {

        final ConcurrentLinkedQueue<SSC> sscCollection = new ConcurrentLinkedQueue<>();
        final List<Callable<List<SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (final DataSet dataSet : dataSets) {
            callables.add(() -> Fragmentation.buildSSCCollection(dataSet.getStructure()
                                                                        .toAtomContainer(), dataSet.getSpectrum(),
                                                                 dataSet.getAssignment(), maxSphere));
        }
        ParallelTasks.processTasks(callables, sscCollectionTemp -> sscCollectionTemp.parallelStream()
                                                                                    .forEach(sscCollection::add),
                                   nThreads);

        return new ArrayList<>(sscCollection);
    }

    /**
     * Builds a collection of substructure-subspectrum-correlations (SSC objects) from one
     * structure for all its atoms by using a breadth first search
     * with spherical limit.
     *
     * @param structure  Structure to build the SSCs from
     * @param spectrum   Spectrum with signals to use for each assigned atom in structure
     * @param assignment Assignments between atoms in structure and belonging signals in spectrum
     * @param maxSphere  Spherical limit for building a substructure into
     *                   all directions to be the same type as in used spectrum property.
     *
     * @return
     *
     * @see Fragmentation#buildSSC(IAtomContainer, Spectrum, Assignment, int, int)
     */
    public static List<SSC> buildSSCCollection(final IAtomContainer structure, final Spectrum spectrum,
                                               final Assignment assignment, final int maxSphere) throws CDKException {
        // if the structure contains explicit hydrogens atoms then, after
        // removing them, the assignments have to be corrected
        if (casekit.nmr.Utils.containsExplicitHydrogens(structure)) {
            //            final Map<IAtom, Integer> prevAtomIndices = casekit.nmr.Utils.convertExplicitToImplicitHydrogens(structure);
            //            for (final IAtom atom : prevAtomIndices.keySet()) {
            //                // consider only atoms with same atom type as first spectrum nucleus type
            //                if (!atom.getSymbol()
            //                         .equals(Utils.getAtomTypeFromSpectrum(spectrum, 0))) {
            //                    continue;
            //                }
            //                for (int i = 0; i
            //                        < structure.getAtomCount(); i++) {
            //                    if ((structure.getAtom(i)
            //                            == atom)
            //                            && (i
            //                            != prevAtomIndices.get(atom))) {
            //                        assignment.setAssignment(0, assignment.getIndex(0, prevAtomIndices.get(atom)), new int[]{i});
            //                        break;
            //                    }
            //                }
            //            }

            return new ArrayList<>();
        }
        final List<SSC> sscCollection = new ArrayList<>();
        SSC ssc;
        for (int i = 0; i
                < structure.getAtomCount(); i++) {
            ssc = Fragmentation.buildSSC(structure, spectrum, assignment, i, maxSphere);
            // if one part of the structure could not be built then skip the whole structure
            // and return an empty SSC library
            if (ssc
                    == null) {
                return new ArrayList<>();
            }
            sscCollection.add(ssc);
        }

        return sscCollection;
    }

    /**
     * Builds a substructure-subspectrum-correlation ({@link model.SSC}) object
     * from a structure, a spectrum and signal to atom assignments.
     * The structure fragmentation is done by using breadth first
     * search with a spherical limit and each atom as starting point.
     *
     * @param structure     structure to fragment into substructures
     * @param spectrum      spectrum to split into subspectra
     * @param assignment    signal to atom assignments
     * @param rootAtomIndex Index of start atom
     * @param maxSphere     Spherical limit for building a substructure into
     *                      all directions
     *
     * @return
     *
     * @throws org.openscience.cdk.exception.CDKException
     * @see Fragmentation#buildSubstructure(org.openscience.cdk.interfaces.IAtomContainer, int, int)
     */
    public static SSC buildSSC(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment,
                               final int rootAtomIndex, final int maxSphere) throws CDKException {
        casekit.nmr.Utils.setAromaticityAndKekulize(structure);
        final List<Integer> substructureAtomIndices = Fragmentation.buildSubstructureAtomIndicesList(structure,
                                                                                                     rootAtomIndex,
                                                                                                     maxSphere);
        final IAtomContainer substructure = Fragmentation.buildSubstructure(structure, rootAtomIndex, maxSphere);

        final Spectrum subspectrum = new Spectrum();
        subspectrum.setNuclei(spectrum.getNuclei());
        subspectrum.setSignals(new ArrayList<>());
        final Assignment subassignment = new Assignment();
        subassignment.setNuclei(spectrum.getNuclei());
        subassignment.initAssignments(subspectrum.getSignalCount());
        IAtom atomInStructure;
        Signal signalToAdd;
        for (int j = 0; j
                < substructure.getAtomCount(); j++) {
            atomInStructure = structure.getAtom(substructureAtomIndices.get(j));
            if (atomInStructure.getSymbol()
                               .equals(Utils.getAtomTypeFromSpectrum(subspectrum, 0))) {
                if ((assignment.getIndex(0, substructureAtomIndices.get(j))
                        == null)
                        || (spectrum.getSignal(assignment.getIndex(0, substructureAtomIndices.get(j)))
                        == null)) {
                    return null;
                }
                signalToAdd = spectrum.getSignal(assignment.getIndex(0, substructureAtomIndices.get(j)))
                                      .buildClone();
                AssemblyUtils.addSignalToSSC(subspectrum, subassignment, signalToAdd, j);
            }
        }
        subspectrum.setSolvent(spectrum.getSolvent());
        subspectrum.setSpectrometerFrequency(spectrum.getSpectrometerFrequency());

        return new SSC(substructure, subspectrum, subassignment, 0, utils.Utils.buildHOSECodes(substructure, maxSphere),
                       utils.Utils.getUnsaturatedAtomIndices(substructure));
    }

    /**
     * Builds a substructure from a structure using a breadth first search
     * with spherical limit, starting point as well as HOSE code priority order
     * of next neighbor atoms.
     *
     * @param structure     IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxSphere     Spherical limit for building a substructure into
     *                      all directions
     *
     * @return
     *
     * @see HOSECodeBuilder#buildConnectionTree(IAtomContainer, int, Integer)
     * @see HOSECodeBuilder#buildAtomContainer(ConnectionTree)
     */
    public static IAtomContainer buildSubstructure(final IAtomContainer structure, final int rootAtomIndex,
                                                   final int maxSphere) {
        return HOSECodeBuilder.buildAtomContainer(
                HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxSphere));
    }

    /**
     * Builds a set of substructure atom indices from a structure using a
     * breadth first search with spherical limit, starting point as well as
     * HOSE code priority order of next neighbor atoms.
     *
     * @param structure     IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxSphere     Spherical limit for building a substructure into
     *                      all directions
     *
     * @return
     *
     * @see HOSECodeBuilder#buildConnectionTree(IAtomContainer, int, Integer)
     */
    public static ArrayList<Integer> buildSubstructureAtomIndicesList(final IAtomContainer structure,
                                                                      final int rootAtomIndex, final int maxSphere) {
        return HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxSphere)
                              .getKeys();
    }
}
