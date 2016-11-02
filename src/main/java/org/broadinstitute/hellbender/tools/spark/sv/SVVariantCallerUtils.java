package org.broadinstitute.hellbender.tools.spark.sv;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;
import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion.InversionType.*;

/**
 * Variants utility functions helping calling structural variants.
 */
final class SVVariantCallerUtils {

    static boolean isInversion(final BreakpointAllele allele) {
        try {
            final BreakpointAlleleInversion invAllele = (BreakpointAlleleInversion) allele;
            return invAllele.leftAlignedLeftBreakpoint.getContig().equals(invAllele.leftAlignedRightBreakpoint.getContig())
                    &&
                    (invAllele.getInversionType() == INV_3_TO_5 || invAllele.getInversionType() == INV_5_TO_3);
        } catch (final ClassCastException ccex) {
            return false;
        }
    }
}
