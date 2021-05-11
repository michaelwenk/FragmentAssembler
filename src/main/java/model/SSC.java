package model;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.Spectrum;
import lombok.*;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.List;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
public class SSC {
    private IAtomContainer structure;
    private Spectrum spectrum;
    private Assignment assignment;
    private List<String> hoseCodes;
    private List<Integer> unsaturatedAtomIndices;

    public SSC buildClone() {
        try {
            return new SSC(this.structure.clone(), this.getSpectrum()
                                                       .buildClone(), this.getAssignment()
                                                                          .buildClone(),
                           new ArrayList<>(this.getHoseCodes()), new ArrayList<>(this.getUnsaturatedAtomIndices()));
        } catch (final CloneNotSupportedException e) {
            e.printStackTrace();
        }

        return null;
    }
}
