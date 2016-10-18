package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public final class ReadsContextUnitTest extends BaseTest {

    @DataProvider(name = "EmptyReadsContextDataProvider")
    public Object[][] getEmptyReadsContextData() {
        // Default-constructed ReadsContexts and ReadsContexts constructed from null ReadsDataSources/intervals
        // should behave as empty context objects.
        return new Object[][] {
                { new ReadsContext() },
                { new ReadsContext(null, null) },
                { new ReadsContext(null, new SimpleInterval("1", 1, 1) ) },
                { new ReadsContext(new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")), null) }
        };
    }

    @Test(dataProvider = "EmptyReadsContextDataProvider")
    public void testEmptyReadsContext( final ReadsContext readsContext ) {
        Assert.assertFalse(readsContext.hasBackingDataSource() && readsContext.getInterval() != null,
                           "Empty ReadsContext reports having both a backing data source and an interval");
        Assert.assertFalse(readsContext.iterator().hasNext(), "Empty ReadsContext should have returned an empty iterator from iterator()");
    }

    @DataProvider(name = "ValidReadsContextDataProvider")
    public Object[][] getValidReadsContextData() {
        // Default-constructed ReadsContexts and ReadsContexts constructed from null ReadsDataSources/intervals
        // should behave as empty context objects.
        return new Object[][]{
                {new ReadsContext(
                        new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")),
                        new SimpleInterval("1", 200, 210),
                        ReadFilterLibrary.ALLOW_ALL_READS), true
                },
                {new ReadsContext(
                        new ReadsDataSource(new File(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam")),
                        new SimpleInterval("1", 200, 210),
                        new ReadFilter() {
                            private static final long serialVersionUID = 1L;
                            @Override public boolean test(final GATKRead read) {return false;}}), false
                }
        };
    }

    @Test(dataProvider = "ValidReadsContextDataProvider")
    public void testValidReadsContext(final ReadsContext readsContext, final boolean hasReads ) {
        Assert.assertTrue(readsContext.hasBackingDataSource() && readsContext.getInterval() != null,
                "Valid ReadsContext reports having no backing data source or interval");
        Assert.assertEquals(readsContext.iterator().hasNext(), hasReads, "Valid ReadsContext should have returned a non empty iterator from iterator()");
    }

}