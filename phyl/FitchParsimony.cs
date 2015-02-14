using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


    /**
    * Class FitchParsimony provides the Fitch algorithm for computing the parsimony
    * score of a {@linkplain DnaSequenceTree}. For further information, see:
    * <UL>
    * <LI>
    * W. Fitch. Toward defining the course of evolution: minimum change for a
    * specified tree topology. <I>Systematic Zoology,</I> 20:406-416, 1971.
    * <LI>
    * J. Felsenstein. <I>Inferring Phylogenies.</I> Sinauer Associates, 2004, pages
    * 11-13.
    * </UL>
    *
    * @author Clark Wang 
    * @version 18-Feb-2015
    */


namespace phyl
{
    class FitchParsimony
    {

        // Prevent construction.
        private FitchParsimony()
        {
        }

        public static int computeScore (DnaSequenceTree tree)
        {
            int root = tree.root();
            computeScore(tree, root);
            return tree.seq(root).score();
        }

        /**
        * Compute the Fitch parsimony score of the given node in the given DNA
        * sequence tree.
        *
        * @param tree DNA sequence tree.
        * @param index Node index.
        */
        private static void computeScore (DnaSequenceTree tree,int index)
        {
            // Stop recursion at a tip node.
            int child1 = tree.child1(index);
            int child2 = tree.child2(index);
            if (child1 == -1) return;
            // Compute scores of child nodes.
            computeScore(tree, child1);
            computeScore(tree, child2);
            // Associate a new DNA sequence with this node if necessary.
            DnaSequence seq1 = tree.seq(child1);
            DnaSequence seq2 = tree.seq(child2);
            DnaSequence seq = tree.seq(index);
            if (seq == null)
            {
                seq = new DnaSequence(seq1.length());
                tree.seq(index, seq);
            }
            // Set this node's DNA sequence to the Fitch ancestor of the two child
            // nodes' DNA sequences.
            seq.setFitchAncestor(seq1, seq2);
            seq.name("" + (seq.score() - seq1.score() - seq2.score()));
        }


        public static int updateScore (DnaSequenceTree tree,int tip,DnaSequence[] seqarray)
        {
            int i = 0;
            // Update all nodes from tip's parent through root.
            DnaSequence seq = tree.seq(tip);
            int index = tree.parent(tip);
            while (index != -1)
            {
                int child1 = tree.child1(index);
                int child2 = tree.child2(index);
                DnaSequence seq1 = tree.seq(child1);
                DnaSequence seq2 = tree.seq(child2);
                seq = seqarray[i++];
                seq.setFitchAncestor(seq1, seq2);
                tree.seq(index, seq);
                index = tree.parent(index);
            }
            return seq.score();
        }

    }
}
