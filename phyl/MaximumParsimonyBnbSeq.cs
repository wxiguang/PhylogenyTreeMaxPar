using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phyl
{
    class MaximumParsimonyBnbSeq
    {

        // List of DNA sequences with which to construct trees.
        private DnaSequenceList seqList;
        // Initial bound.
        private int initialBound;
        // For holding search results.
        private MaximumParsimonyResults results;
        // Length of each DNA sequence.
        private int L;
        // Number of DNA sequences.
        private int N;
        // Tree capacity.
        private int C;
        // Number of absent states as each DNA sequence is added.
        private int[] absentStates;
        // Stack of DNA sequence trees.
        private DnaSequenceTree[] treeStack;
        // Stack of auxiliary DNA sequence arrays.
        DnaSequence[][] seqArrayStack;
        // Tree signature currently being searched.
        private int[] signature;



        /**
        * Construct a new maximum parsimony phylogenetic tree construction
        * algorithm object.
        *
        * @param seqList DNA sequence list.
        * @param initialBound Initial bound for branch-and-bound search.
        * @param results Object in which to store the results.
        */

        	public MaximumParsimonyBnbSeq (DnaSequenceList seqList,int initialBound,MaximumParsimonyResults results)
            {
                // Record parameters.
                this.seqList = seqList;
                this.initialBound = initialBound;
                this.results = results;
                // Initialize.
                L = seqList.seq(0).length();
                N = seqList.length();
                C = 2*N - 1;
                // Compute number of absent states as each DNA sequence is added.
                absentStates = seqList.countAbsentStates();
                // Set up stack of DNA sequence trees.
                treeStack = new DnaSequenceTree [N];
                for (int i = 0; i < N; ++ i)
                {
                treeStack[i] = new DnaSequenceTree (C);
                }
                // Initialize DNA sequence tree at first level of the search graph.
                treeStack[0].add(0, seqList.seq(0));
                // Set up stack of auxiliary DNA sequence arrays.
                seqArrayStack = new DnaSequence[N][];
                for (int i = 0; i < N; ++i)
                {
                    DnaSequence[] seqArray = new DnaSequence[i];
                    seqArrayStack[i] = seqArray;
                    for (int j = 0; j < i; ++j)
                    {
                        seqArray[j] = new DnaSequence(L);
                    }
                }
                // Set up tree signature.
                signature = new int[N];
            }
        /**
        * Find the maximum parsimony phylogenetic tree(s) in the search graph. The
        * DNA sequence list was specified to the constructor. The results are
        * stored in the {@linkplain MaximumParsimonyResults} object specified to
        * the constructor. The <TT>findTrees()</TT> method will only find trees
        * whose parsimony scores are less than or equal to the
        * <TT>initialBound</TT> specified to the constructor or the best bound
        * found thereafter, whichever is smaller.
        */
        	public void findTrees() 
            {
                // Initialize tree signature.
                signature[0] = 0;
                for (int i = 1; i < N; ++ i)
                {
                signature[i] = -1;
                }
                // Traverse remaining levels of the search graph.
                int level = 1;
                bool done = false;
                results.clear();
                results.score =initialBound;
                while (! done)
                {
                DnaSequenceTree prevTree = treeStack[level-1];
                // If we have reached the bottom of the search graph, we have a
                // tentative solution.
                if (level == N)
                {
                int tentativeScore = prevTree.seq (prevTree.root()) .score();
                // Record tentative solution.
                results.add (signature, tentativeScore);
                // Go to previous level.
                -- level;
                if (level == 1) done = true;
                }
                // If there are no more positions to try at this level, reset
                // position at this level and go to previous level.
                else if (signature[level] == 2 * (level - 1))
                {
                    signature[level] = -1;
                    --level;
                    if (level == 1) done = true;
                }
                // If there are more positions to try at this level, add the DNA
                // sequence to the tree at the next position and do
                // branch-and-bound.
                else
                {
                    ++signature[level];
                    DnaSequenceTree currTree = treeStack[level];
                    currTree.copy(prevTree);
                    int tip = currTree.add(signature[level], seqList.seq(level));
                    int partialScore = FitchParsimony.updateScore(currTree, tip, seqArrayStack[level]);
                    // If partial parsimony score plus number of absent states in
                    // the remaining levels is less than or equal to the best
                    // solution's score, go to the next level (branch), otherwise
                    // try the next choice at this level (bound).
                    if (partialScore + absentStates[level] <= results.score)
                    {
                        ++level;
                    }
                }
                }
            }


    }
}
