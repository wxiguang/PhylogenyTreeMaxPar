using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phyl
{
/**
 * Class JukesCantorDistance provides an object that computes the Jukes-Cantor
 * distance between two {@linkplain DnaSequence}s. This is the corrected
 * distance under the Jukes-Cantor model of DNA sequence evolution. The formula
 * is
 * <CENTER>
 * <I>D</I><SUB><I>JC</I></SUB> = &minus;3/4 <I>N</I> ln (1 &minus; 4/3 <I>D</I><SUB><I>H</I></SUB>/<I>N</I>)
 * </CENTER>
 * where <I>D</I><SUB><I>JC</I></SUB> is the Jukes-Cantor distance,
 * <I>D</I><SUB><I>H</I></SUB> is the Hamming distance (number of differing
 * sites), and <I>N</I> is the number of sites. For further information, see:
 * <UL>
 * <LI>
 * T. Jukes and C. Cantor. Evolution of protein molecules. In M. Munro,
 * editor. <I>Mammalian Protein Metabolism, Volume III.</I> Academic Press,
 * 1969, pages 21-132.
 * <LI>
 * J. Felsenstein. <I>Inferring Phylogenies.</I> Sinauer Associates, 2004,
 * pages 156-158.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 23-Jul-2008
 */
public class JukesCantorDistance:Distance
	{

// Exported constructors.

	/**
	 * Construct a new Jukes-Cantor distance object.
	 */
	public JukesCantorDistance()
		{
		}

// Exported operations.

	/**
	 * Compute the distance between the two given DNA sequences. It is assumed
	 * that the DNA sequences are the same length.
	 *
	 * @param  seq1  First DNA sequence.
	 * @param  seq2  Second DNA sequence.
	 *
	 * @return  Distance.
	 */
	public double distance (DnaSequence seq1,DnaSequence seq2)
		{
		double D = seq1.distance (seq2);
		double N = seq1.length();
		double x = 1.0 - D/N/0.75;
		return x <= 0.0 ? Double.PositiveInfinity :Math.Abs (-0.75*N*Math.Log(x));
		}

	}
}
