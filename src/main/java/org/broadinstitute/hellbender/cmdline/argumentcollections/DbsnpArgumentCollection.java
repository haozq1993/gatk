package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.engine.FeatureInput;

public final class DbsnpArgumentCollection implements ArgumentCollectionDefinition{
    private static final long serialVersionUID = 1L;

    /**
     * A dbSNP VCF file.
     */
    @Argument(fullName="dbsnp", shortName = "D", doc="dbSNP file", optional=true)
    public FeatureInput<VariantContext> dbsnp;

}

