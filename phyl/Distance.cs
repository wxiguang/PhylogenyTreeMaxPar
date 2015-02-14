using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phyl
{
    /**
    * Interface Distance specifies the interface for an object that computes the
    * distance between two {@linkplain DnaSequence}s.
    *
    * @author Alan Kaminsky
    * @version 23-Jul-2008
    */
    public interface Distance
    {
        // Exported operations.
        /**
        * Compute the distance between the two given DNA sequences. It is assumed
        * that the DNA sequences are the same length.
        *
        * @param seq1 First DNA sequence.
        * @param seq2 Second DNA sequence.
        *
        * @return Distance.
        */
        double distance(DnaSequence seq1,DnaSequence seq2);
    }
}
