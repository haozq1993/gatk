package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Unit tests for SVKmer and SVKmerizer.
 */
public class SVKmerUnitTest extends BaseTest {
    @Test
    public void testDefaultConstruction() {
        Assert.assertEquals(new SVKmer(10).toString(10), "AAAAAAAAAA");
        Assert.assertEquals(new SVKmer(11).toString(11), "AAAAAAAAAAA");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooLargeK() {
        final SVKmer tooBigK = new SVKmer(64);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testDefaultConstructionWithTooSmallK() {
        final SVKmer tooSmallK = new SVKmer(0);
    }

    @DataProvider(name = "sequenceStrings")
    public Object[][] getSequenceStrings() {
        return new Object[][] {
                {"ACGTACGTACGT"}, {"ACGTACGTACGTC"}, {"ACGTACGTACGTCC"}, {"ACGTACGTACGTCCC"}, {"ACGTACGTACGTCCCC"}
        };
    }

    @Test(dataProvider = "sequenceStrings")
    public void testConstructionAndToString( final String str ) {
        Assert.assertEquals(str, SVKmerizer.toKmer(str,new SVKmerSmall(31)).toString(str.length()));
        Assert.assertEquals(str, SVKmerizer.toKmer(str.getBytes(),new SVKmerSmall(31)).toString(str.length()));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testSuccessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.delete(0, 1);
        sb.append('A');
        final SVKmer kkk = SVKmerizer.toKmer(str,new SVKmer(str.length()));
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(K-1, 'C');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(K-1, 'G');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(K-1, 'T');
        Assert.assertEquals(sb.toString(), kkk.successor(SVKmer.Base.T, K).toString(K));
    }

    @Test(dataProvider = "sequenceStrings")
    public void testPredecessor( final String str ) {
        final int K = str.length();
        final StringBuilder sb = new StringBuilder(str);
        sb.insert(0, 'A');
        sb.setLength(K);
        final SVKmer kkk = SVKmerizer.toKmer(str,new SVKmer(str.length()));
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.A, K).toString(K));
        sb.setCharAt(0, 'C');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.C, K).toString(K));
        sb.setCharAt(0, 'G');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.G, K).toString(K));
        sb.setCharAt(0, 'T');
        Assert.assertEquals(sb.toString(), kkk.predecessor(SVKmer.Base.T, K).toString(K));
    }

    @Test
    public void testReverseComplementation() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGA",new SVKmer(8)).reverseComplement(8), SVKmerizer.toKmer("TCGTACGT",new SVKmer(8)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmer(9)).reverseComplement(9), SVKmerizer.toKmer("GACGTACGT",new SVKmer(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTCCGTC",new SVKmer(9)).reverseComplement(9), SVKmerizer.toKmer("GACGGACGT",new SVKmer(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTGCGTC",new SVKmer(9)).reverseComplement(9), SVKmerizer.toKmer("GACGCACGT",new SVKmer(9)));
        Assert.assertEquals(SVKmerizer.toKmer("ACGTTCGTC",new SVKmer(9)).reverseComplement(9), SVKmerizer.toKmer("GACGAACGT",new SVKmer(9)));
    }

    @Test
    public void testCanonicalization() {
        Assert.assertEquals(SVKmerizer.toKmer("ACGTACGTC",new SVKmer(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmer(9)));
        Assert.assertEquals(SVKmerizer.toKmer("GACGTACGT",new SVKmer(9)).canonical(9), SVKmerizer.toKmer("ACGTACGTC",new SVKmer(9)));
    }

    @Test
    public void testComparison() {
        final SVKmer kkk1 = SVKmerizer.toKmer("ACGTA",new SVKmer(4));
        final SVKmer kkk2 = SVKmerizer.toKmer("ACGTC",new SVKmer(4));
        Assert.assertTrue(kkk1.compareTo(kkk1) == 0);
        Assert.assertTrue(kkk2.compareTo(kkk2) == 0);
        Assert.assertTrue(kkk1.compareTo(kkk2) < 0);
        Assert.assertTrue(kkk2.compareTo(kkk1) > 0);
    }

    @Test
    public void testHashCode() {
        Assert.assertNotEquals(SVKmerizer.toKmer("TAGCGTA",new SVKmer(7)).hashCode(), SVKmerizer.toKmer("TAGCGTC",new SVKmer(7)).hashCode());
        Assert.assertEquals(SVKmerizer.toKmer("TAGGGTC",new SVKmer(7)).hashCode(), SVKmerizer.toKmer("TAGGGTC",new SVKmer(7)).hashCode());
    }

    @Test
    public void testKmerization() {
        final SVKmerizer<SVKmer> kmerizer = new SVKmerizer<SVKmer>("AAAAATT", 5, new SVKmer(7));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmer(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAT",new SVKmer(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAATT",new SVKmer(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }

    @Test
    public void testKmerizationAcrossN() {
        final SVKmerizer<SVKmer> kmerizer = new SVKmerizer<SVKmer>("AAAAANTTTTT", 5, new SVKmer(11));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("AAAAA",new SVKmer(5)));
        Assert.assertTrue(kmerizer.hasNext());
        Assert.assertEquals(kmerizer.next(), SVKmerizer.toKmer("TTTTT",new SVKmer(5)));
        Assert.assertTrue(!kmerizer.hasNext());
    }
}
