using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

/**
* Class DnaSequence encapsulates a DNA sequence. The DNA sequence consists of a
* sequence of <B>sites</B>. Each site has a <B>state,</B> which is a set of
* <B>bases</B>. The four bases are adenine, cytosine, guanine, and thymine. For
* textual I/O, each state is represented by a single character as follows:
*  * The DNA sequence has an associated <B>score,</B> an integer. The score can be
* set to anything and later retrieved.
* <P>
* The DNA sequence has a <B>name,</B> a string. The name can be set to anything
* and later retrieved.
*
* @author Clark Wang
* @version 10-Feb-2015
*/

namespace phyl
{
    public class DnaSequence
    {
        // Amount of extra padding in byte array.
        private static int PAD = 128;
        // Mapping from the state of a site to the corresponding output character.
        // A=1, C=2, G=4, T=8.
        	static char[] state2char = new char[]
            {/*----*/ '-',
            /*---A*/ 'A',
            /*--C-*/ 'C',
            /*--CA*/ 'M',
            /*-G--*/ 'G',
            /*-G-A*/ 'R',
            /*-GC-*/ 'S',
            /*-GCA*/ 'V',
            /*T---*/ 'T',
            /*T--A*/ 'W',
            /*T-C-*/ 'Y',
            /*T-CA*/ 'H',
            /*TG--*/ 'K',
            /*TG-A*/ 'D',
            /*TGC-*/ 'B',
            /*TGCA*/ 'X'};
        	// Mapping from the state of a site to the number of bits turned on at that
            // site. A=1, C=2, G=4, T=8.
        	public static int[] state2bitCount = new int[]
            {/*----*/ 0,
            /*---A*/ 1,
            /*--C-*/ 1,
            /*--CA*/ 2,
            /*-G--*/ 1,
            /*-G-A*/ 2,
            /*-GC-*/ 2,
            /*-GCA*/ 3,
            /*T---*/ 1,
            /*T--A*/ 2,
            /*T-C-*/ 2,
            /*T-CA*/ 3,
            /*TG--*/ 2,
            /*TG-A*/ 3,
            /*TGC-*/ 3,
            /*TGCA*/ 4};

            // Sequence data. Each site's set of bases is stored as a bitmap in one
            // byte. A=1, C=2, G=4, T=8. PAD bytes of padding are added to avert cache
            // interference.
            public byte[] mySites;
            // Length.
            public int myLength;
            // Score.
            public int myScore;
            // Name.
            public string myName;

        	// 128 bytes of extra padding to avert cache interference.
            private long p0, p1, p2, p3, p4, p5, p6, p7;
            private long p8, p9, pa, pb, pc, pd, pe, pf;
            //private NonSerializedAttribute long p8, p9, pa, pb, pc, pd, pe, pf;
            
            /**
            * Construct a new zero-length DNA sequence. The score is initially 0. The
            * name is initially null.
            */
            public DnaSequence()
            {
                this.myLength = 0;
                this.myScore = 0;
                this.myName = null;
            }
            public DnaSequence (int N)
            {
                this.myLength = N;
                this.myScore = 0;
                this.myName = null;
            }
            public DnaSequence(int N,int score)
            {
                this.myLength = N;
                this.myScore = score;
                this.myName = null;
            }
            /**
            * Construct a new DNA sequence with the given length, score, and name.
            *
            * @param N Length (number of sites).
            * @param score Score.
            * @param name Name. May be null.
            *
            * @exception IllegalArgumentException
            * (unchecked exception) Thrown if <TT>N</TT> &lt; 0.
            */
            public DnaSequence (int N, int score, string name)
            {
                if (N < 0)
                {
                    throw new ArgumentException ("DnaSequence(): N (= " + N + ") < 0, illegal)");
                }
                this.mySites = new byte[N + PAD];
                this.myLength = N;
                this.myScore = score;
                this.myName = name;
            }
            /**
            * Construct a new DNA sequence that is a copy of the given DNA sequence.
            *
            * @param seq DNA sequence to copy.
            *
            * @exception NullPointerException
            * (unchecked exception) Thrown if <TT>seq</TT> is null.
            */
            public DnaSequence(DnaSequence seq)
            {
                this.mySites = seq.mySites;
                this.myLength = seq.myLength;
                this.myScore = seq.myScore;
                this.myName = seq.myName;
            }
            /**
            * Get this DNA sequence's length.
            *
            * @return Length (number of sites).
            */
            public int length()
            {
                return myLength;
            }
            /**
            * Get this DNA sequence's score.
            *
            * @return Score.
            */
            public int score()
            {
                return myScore;
            }
            /**
            * Set this DNA sequence's score.
            *
            * @param score Score.
            */
            public void score (int score)
            {
                myScore = score;
            }
            /**
            * Get this DNA sequence's name.
            *
            * @return Name. May be null.
            */
            public string name()
            {
                return myName;
            }
            /**
            * Set this DNA sequence's name.
            *
            * @param name Name. May be null.
            */
            public void name (string name)
            {
                myName = name;
            }
            /**
            * Make this DNA sequence's sites be the same as the given DNA sequence. It
            * is assumed that this DNA sequence and the given DNA sequence are the same
            * length. This DNA sequence's score and name are unchanged.
            *
            * @param seq DNA sequence to copy.
            *
            * @exception NullPointerException
            * (unchecked exception) Thrown if <TT>seq</TT> is null.
            */
            public void copySites (DnaSequence seq)
            {
                Array.Copy(seq.mySites, 0, this.mySites, 0, myLength);
            }
            /// <summary>
            /// ////////////////////////////////////////////////////////////
            /// </summary>
            /// <param name="o"></param>
            /// <returns></returns>
            public bool equals(Object o)
            {
                if (this == o) return true;
                if ((o == null) || (o.GetType() != this.GetType())) return false;
                DnaSequence ds = (DnaSequence)o;
                if (myLength != ds.myLength) return false;
                if (myScore != ds.myScore) return false;
                if (!(((myName != null) && myName.Equals(ds.myName))
                || ((myName == null) && (ds.myName == null))))
                    return false;
                for (int i = 0; i < myLength; i++)
                {
                    if (mySites[i] != ds.mySites[i])
                        return false;
                }
                return true;
            }

            public int hashCode() {
                int ret = 1669;
                ret = 709 * ret + myLength;
                ret = 709 * ret + myScore;
                ret = 709 * ret + mySites.GetHashCode();
                if (myName != null) {
                }
                return ret;
            }


            /**
            * Compute the distance between this DNA sequence and the given DNA
            * sequence. It is assumed that this DNA sequence and the given DNA sequence
            * are the same length. The distance is the number of differing sites
            * between the two sequences (the Hamming distance).
            *
            * @param seq DNA sequence.
            *
            * @return Distance.
            */
            public double distance (DnaSequence seq)
            {
                byte[] site1 = this.mySites;
                byte[] site2 = seq.mySites;
                int diff = 0;
                int N = myLength;
                for (int i = 0; i < N; ++i)
                {
                    if (site1[i] != site2[i]) ++diff;
                }
                return diff;
            }

            /**
            * Make this DNA sequence be the ancestor of the two given DNA sequences in
            * the Fitch parsimony score algorithm. This DNA sequence's sites are set
            * based on <TT>seq1</TT>'s and <TT>seq2</TT>'s sites. It is assumed that
            * this DNA sequence and the given DNA sequences are the same length. This
            * DNA sequence's score is set to the sum of <TT>seq1</TT>'s score,
            * <TT>seq2</TT>'s score, and the number of state changes at the ancestor.
            * This DNA sequence's name is unchanged.
            *
            * @param seq1 First child DNA sequence.
            * @param seq2 Second child DNA sequence.
            */
            public void setFitchAncestor (DnaSequence seq1,DnaSequence seq2)
            {
                // Get references to sites.
                byte[] ancestor = this.mySites;
                byte[] descendent1 = seq1.mySites;
                byte[] descendent2 = seq2.mySites;
                int N = myLength;
                // Process all sites. Count state changes.
                int nChanges = 0;
                for (int i = 0; i < N; ++i)
                {
                    // Compute intersection of states.
                    int state1 = descendent1[i];
                    int state2 = descendent2[i];
                    int state3 = state1 & state2;
                    // If intersection is not empty, record intersection, otherwise
                    // record union and note one state change.
                    if (state3 == 0)
                    {
                        state3 = state1 | state2;
                        ++nChanges;
                    }
                    // Update site.
                    ancestor[i] = (byte)state3;
                }
                // Record number of state changes.
                this.myScore = seq1.myScore + seq2.myScore + nChanges;
            }
            /**
            * Returns a string version of this DNA sequence. The string consists of
            * just the sequence of states (the score and name are not included).
            */
            public string toString()
            {
                StringBuilder buf = new StringBuilder();
                byte[] site = mySites;
                int N = myLength;
                for (int i = 0; i < N; ++i)
                {
                    buf.Append(state2char[site[i]]);
                }
                return buf.ToString();
            }


    }
}
