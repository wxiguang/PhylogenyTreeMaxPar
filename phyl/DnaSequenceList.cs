using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;


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


    class Scanner : System.IO.StringReader
    {
        string currentWord;

        public Scanner(string source)
            : base(source)
        {
            readNextWord();
        }

        private void readNextWord()
        {
            System.Text.StringBuilder sb = new StringBuilder();
            char nextChar;
            int next;
            do
            {
                next = this.Read();
                if (next < 0)
                    break;
                nextChar = (char)next;
                if (char.IsWhiteSpace(nextChar))
                    break;
                sb.Append(nextChar);
            } while (true);
            while ((this.Peek() >= 0) && (char.IsWhiteSpace((char)this.Peek())))
                this.Read();
            if (sb.Length > 0)
                currentWord = sb.ToString();
            else
                currentWord = null;
        }

        public bool hasNextInt()
        {
            if (currentWord == null)
                return false;
            int dummy;
            return int.TryParse(currentWord, out dummy);
        }

        public int nextInt()
        {
            try
            {
                return int.Parse(currentWord);
            }
            finally
            {
                readNextWord();
            }
        }

        public bool hasNextDouble()
        {
            if (currentWord == null)
                return false;
            double dummy;
            return double.TryParse(currentWord, out dummy);
        }

        public double nextDouble()
        {
            try
            {
                return double.Parse(currentWord);
            }
            finally
            {
                readNextWord();
            }
        }

        public bool hasNext()
        {
            return currentWord != null;
        }
    }
    
    
    public static class ExtensionMethod
    {
        public static void ArraysFill<T>(this T[] originalArray, T with)
        {
            for (int i = 0; i < originalArray.Length; i++)
            {
                originalArray[i] = with;
            }
        }
    }
    [Serializable]
    public class DnaSequenceList 
    {
        // DNA sequences.
        public DnaSequence[] mySequence;
        // Mapping from site (index) to whether site is informative (true/false). If
        // null, must be recomputed.
        private bool[] isInformative;
        // Number of informative sites.
        private int nInformative;
        // Number of state changes in uninformative sites.
        private int nChanges;

        /**
        * Construct a new DNA sequence list.
        */
        public DnaSequenceList()
        {
        }
        /**
        * Construct a new DNA sequence list that is a copy of the given DNA
        * sequence list.
        * <P>
        * <I>Note:</I> The DNA sequences in the new list are copies of (not
        * references to) the DNA sequences in the given list.
        *
        * @param list DNA sequence list to copy.
        *
        * @exception NullPointerException
        * (unchecked exception) Thrown if <TT>list</TT> is null.
        */
        public DnaSequenceList(DnaSequenceList list)
        {
            int N = list.mySequence.Length;
            this.mySequence = new DnaSequence[N];
            for (int i = 0; i < N; ++i)
            {
                this.mySequence[i] = new DnaSequence(list.mySequence[i]);
            }
            if (list.isInformative != null)
            {
                this.isInformative = (bool[])list.isInformative.Clone();
            }
            this.nInformative = list.nInformative;
            this.nChanges = list.nChanges;
        }
        /**
        * Obtain this DNA sequence list's length.
        *
        * @return Length <I>N</I> (number of DNA sequences).
        */
        public int length()
        {
            return mySequence.Length;
        }

        /**
        * Get the DNA sequence at the given index in this DNA sequence list.
        *
        * @param i Index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
        *
        * @return DNA sequence.
        *
        * @exception ArrayIndexOutOfBoundsException
        * (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
        */
        public DnaSequence seq(int i)
        {
            return mySequence[i];
        }


        	/**
	 * Read a DNA sequence list from the given input file. The input file must
	 * be in interleaved PHYLIP format.
	 * <P>
	 * The DNA sequences' sites and names are read from the input file. The DNA
	 * sequences' scores are set to 0.
	 *
	 * @param  file  File.
	 *
	 * @return  DNA sequence list.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred. Thrown if the input file's contents
	 *     were invalid.
	 */
        public static DnaSequenceList read()
        {

            string path = @"C:\Users\Administrator.USER-20141204TC.000\Desktop\dragon.phy";
            string file = "dragon.phy";
            //Scanner filescanner = new Scanner(file);
            FileStream filescanner = new FileStream(path, FileMode.Open, FileAccess.Read);
            //Scanner linescanner;
            StreamReader sr = new StreamReader(filescanner);
            //使用StreamReader类来读取文件 
            sr.BaseStream.Seek(0, SeekOrigin.Begin);
            Scanner linescanner = new Scanner(path);
            int S, N;
            DnaSequenceList list;
            int[] sitecount;
            string line;
            
            try
            {
                
                // Read number of species and number of sites from first line.
                if (!filescanner.CanRead)
                {
                    throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Empty file");
                }
                linescanner = new Scanner(sr.ReadLine());
                if (!linescanner.hasNextInt())
                {
                    throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Number of species invalid or missing");
                }
                S = linescanner.nextInt();
                if (S < 2)
                {
                    throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Number of species must be >= 2");
                }
                if (!linescanner.hasNextInt())
                {
                    throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Number of sites invalid or missing");
                }
                N = linescanner.nextInt();
                if (N < 1)
                {
                    throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Number of sites must be >= 1");
                }
                
             

                // Set up DNA sequence list and site count array.
                list = new DnaSequenceList();
                list.mySequence = new DnaSequence[S];
                sitecount = new int[S];

                // Read sequence data from groups of S lines until EOF.
                for (; ; )
                {
                    for (int s = 0; s < S; ++s)
                    {
                        // Get a line of sequence data for species s.
                        if (sr.Peek()>=0)
                        {
                        }
                        else if (s != 0 || sitecount[s] == 0)
                        {
                            throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Missing a line of sequence data for species " + (s + 1));
                        }
                        else
                        {
                            goto fileloopBreak;
                        }
                        line = sr.ReadLine();

                        // Ignore blank lines.
                        if (line.Trim().Equals(""))
                        {
                            --s;
                            continue;
                        }

                        // The first time, extract sequence name and create
                        // DnaSequence object.
                        if (sitecount[s] == 0)
                        {
                            if (line.Length < 10)
                            {
                                throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Name must be 10 characters for species " + (s + 1));
                            }
                            list.mySequence[s] = new DnaSequence(N, 0, line.Substring(0, 10).Trim());
                            line = line.Substring(10);
                        }

                        // Parse characters in sequence data.
                        int len = line.Length;
                        byte[] seq_Renamed = list.mySequence[s].mySites;
                        byte[] seq0 = list.mySequence[0].mySites;
                        int count = sitecount[s];
                        for (int i = 0; i < len; ++i)
                        {
                            switch (line[i])
                            {
                                case 'O':
                                case 'o':
                                case '-':
                                    verifyCount(count, N,s);
                                    seq_Renamed[count] = (byte)0; // ----
                                    ++count;
                                    break;
                                case 'A':
                                case 'a':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)1; // ---A
                                    ++count;
                                    break;
                                case 'C':
                                case 'c':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)2; // --C-
                                    ++count;
                                    break;
                                case 'M':
                                case 'm':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)3; // --CA
                                    ++count;
                                    break;
                                case 'G':
                                case 'g':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)4; // -G--
                                    ++count;
                                    break;
                                case 'R':
                                case 'r':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)5; // -G-A
                                    ++count;
                                    break;
                                case 'S':
                                case 's':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)6; // -GC-
                                    ++count;
                                    break;
                                case 'V':
                                case 'v':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)7; // -GCA
                                    ++count;
                                    break;
                                case 'T':
                                case 't':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)8; // T---
                                    ++count;
                                    break;
                                case 'W':
                                case 'w':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)9; // T--A
                                    ++count;
                                    break;
                                case 'Y':
                                case 'y':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)10; // T-C-
                                    ++count;
                                    break;
                                case 'H':
                                case 'h':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)11; // T-CA
                                    ++count;
                                    break;
                                case 'K':
                                case 'k':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)12; // TG--
                                    ++count;
                                    break;
                                case 'D':
                                case 'd':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)13; // TG-A
                                    ++count;
                                    break;
                                case 'B':
                                case 'b':
                                    verifyCount(count, N,s);
                                    seq_Renamed[count] = (byte)14; // TGC-
                                    ++count;
                                    break;
                                case 'X':
                                case 'x':
                                case 'N':
                                case 'n':
                                case '?':
                                    verifyCount(count, N, s);
                                    seq_Renamed[count] = (byte)15; // TGCA
                                    ++count;
                                    break;
                                case '.':
                                    verifyCount(count, N, s);
                                    if (s == 0)
                                    {
                                        throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "'.' not allowed in species 1");
                                    }
                                    if (count >= sitecount[0])
                                    {
                                        throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "'.' in species " + (s + 1) + " has no corresponding site in species 1");
                                    }
                                    seq_Renamed[count] = seq0[count];
                                    ++count;
                                    break;
                            }
                        }
                        sitecount[s] = count;
                    speciesloopContinue: ;
                    }
                speciesloopBreak: ;
                fileloopContinue: ;
                }
            fileloopBreak:

                // Verify correct site count for all species.
                for (int s = 0; s < S; ++s)
                {
                    if (sitecount[s] < N)
                    {
                        throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Too few sites for species " + (s + 1));
                    }
                    else if (sitecount[s] > N)
                    {
                        throw new IOException("DnaSequenceList.read(\"" + file + "\"): " + "Too many sites for species " + (s + 1));
                    }
                }

                // Return DNA sequence list.
                return list;
            }

            finally
            {
                //filescanner.close();
                linescanner.Close();
                filescanner.Close();
            }
        }


        private static void verifyCount(int count, int N, int s)
        {
            if (count >= N)
            {
                throw new IOException("DnaSequenceList.read(\""  + "\"): " + "Too many sites for species " + (s + 1));
            }
        }


        public int exciseUninformativeSites()
        {
            int S = mySequence.Length;
            int N = mySequence[0].length();
            // Determine which sites are informative.
            computeInformativeSites();
            // Excise uninformative sites from sequences.
            for (int s = 0; s < S; ++s)
            {
                byte[] oldSites = mySequence[s].mySites;
                mySequence[s] = new DnaSequence (nInformative,mySequence[s].myScore,mySequence[s].myName);
                byte[] excSites = mySequence[s].mySites;
                int j = 0;
                for (int i = 0; i < N; ++i)
                {
                    if (isInformative[i])
                    {
                        excSites[j++] = oldSites[i];
                    }
                }
            }
            // Mark all sites as informative.
            isInformative = new bool[nInformative];

            for (int i = 0; i < isInformative.Length; i++)
            {
                isInformative[i] = true;
            }
            //ArraysFill(isInformative, true);
            // Return number of state changes.
            return nChanges;
        }

        /**
        * Returns the number of informative sites in this DNA sequence list.
        *
        * @return Number of informative sites.
        */
        public int informativeSiteCount()
        {
            computeInformativeSites();
            return nInformative;
        }

        /**
        * Compute information about informative sites.
        */
        private void computeInformativeSites()
        {
            if (isInformative != null) return;
            int S = mySequence.Length;
            int N = mySequence[0].length();
            // Allocate storage to remember each site's category: true =
            // informative, false = uninformative. Also count number of informative
            // sites and number of state changes in uninformative sites.
            isInformative = new bool[N];
            nInformative = 0;
            nChanges = 0;
            // Allocate storage to count states at each site.
            int[] stateCount = new int[16];
            // Examine all sites.
            for (int i = 0; i < N; ++i)
            {
                //Arrays.fill(stateCount, 0);
                for (int m = 0; m < stateCount.Length; m++)
                {
                    stateCount[m] = 0;
                }

                // Examine current site in all sequences.
                for (int s = 0; s < S; ++s)
                {
                    ++stateCount[mySequence[s].mySites[i]];
                }
                // Count how many values in stateCount are 2 or greater.
                int x = 0;
                for (int j = 0; j < 16; ++j)
                {
                    if (stateCount[j] >= 2) ++x;
                }
                // Categorize current site.
                if (x >= 2)
                {
                    // Informative site.
                    isInformative[i] = true;
                    ++nInformative;
                }
                else
                {
                    // Uninformative site. Increase number of state changes by
                    // (number of different states - 1).
                    isInformative[i] = false;
                    for (int j = 0; j < 16; ++j)
                    {
                        if (stateCount[j] > 0) ++nChanges;
                    }
                    --nChanges;
                }
            }
        }



        /**
         * Determine the number of absent states after adding each sequence in this
         * DNA sequence list to a tree. The return value <I>A</I> is an
         * <I>N</I>-element array, where <I>N</I> is the length of this DNA sequence
         * list. As sequences from this list are added to a tree in order from
         * <I>i</I> = 0 to <I>N</I>&minus;1, <I>A</I>[<I>i</I>] is the number of
         * character states that do not yet appear in the tree. Thus, the number of
         * state changes in the tree must increase by at least <I>A</I>[<I>i</I>]
         * when the sequences after sequence <I>i</I> are added to the tree. This
         * can be used to prune a branch-and-bound search.
         *
         * @return  Array <I>A</I>.
         */
        public int[] countAbsentStates()
        {
            int N = mySequence.Length;
            int L = mySequence[0].length();
            int[] A = new int[N];

            // Compute the union of all the DNA sequences.
            byte[] sites = new byte[L];
            for (int i = 0; i < N; ++i)
            {
                byte[] mysites_i = mySequence[i].mySites;
                for (int j = 0; j < L; ++j)
                {
                    sites[j] |= mysites_i[j];
                }
            }

            // Subtract each sequence from the union, count and record states.
            for (int i = 0; i < N; ++i)
            {
                byte[] mysites_i = mySequence[i].mySites;
                int count = 0;
                for (int j = 0; j < L; ++j)
                {
                    sites[j] &= mysites_i[j];
                    count += DnaSequence.state2bitCount[sites[j]];
                }
                A[i] = count;
            }

            return A;
        }

        /**
         * Create a DNA sequence tree from this DNA sequence list and the given tree
         * signature. The tree signature is an array of indexes of length <I>N</I>,
         * where <I>N</I> is the length of this list. To construct the tree, for all
         * <I>i</I> from 0 to <I>N</I>&minus;1, the DNA sequence at index <I>i</I>
         * in this list is added to the tree at index <TT>signature[i]</TT> using
         * the <TT>DnaSequenceTree.add()</TT> method. For all <I>i</I>,
         * <TT>signature[i]</TT> must be in the range 0 ..
         * 2(<I>i</I>&nbsp;&minus;&nbsp;1), except <TT>signature[0]</TT> is 0.
         * <P>
         * <I>Note:</I> The returned tree has references to (not copies of) the DNA
         * sequences in this list.
         *
         * @param  signature  Tree signature (array of tree indexes).
         *
         * @return  Tree.
         */
        public DnaSequenceTree toTree(int[] signature)
        {
            int N = mySequence.Length;
            DnaSequenceTree tree = new DnaSequenceTree(2 * N - 1);
            for (int i = 0; i < N; ++i)
            {
                tree.add(signature[i], mySequence[i]);
            }
            return tree;
        }




    }





}
