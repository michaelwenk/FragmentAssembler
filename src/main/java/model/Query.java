package model;

import casekit.nmr.model.Spectrum;
import lombok.*;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
public class Query {

    private Spectrum querySpectrum;
    private String mf;
}
