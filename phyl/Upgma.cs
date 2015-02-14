using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phyl
{
    class Upgma
    {
        // Prevent construction.
        private Upgma()
        {
        }

        /**
            * Build a phylogenetic tree of the given DNA sequences. The return value is
            * a phylogenetic tree, including branch lengths, of the DNA sequences in
            * <TT>seqlist</TT>, as constructed by the UPGMA algorithm. The DNA
            * sequences in the tree are references to (not copies of) the DNA sequences
            * in <TT>seqlist</TT>. The distances between the DNA sequences are
            * calculated using the <TT>dcalc</TT> object.
            *
            * @param  seqList  List of DNA sequences.
            * @param  dcalc    Object to calculate distances between DNA sequences.
            *
            * @return  Phylogenetic tree.
            */
        public static DnaSequenceTree buildTree (DnaSequenceList seqList,Distance dcalc)
		{
		// Get initial DNA sequences and put each one in its own tree.
		int N = seqList.length();
		DnaSequenceTree[] tree = new DnaSequenceTree [N+1];
		for (int i = 0; i < N; ++ i)
			{
			tree[i] = new DnaSequenceTree (1);
			tree[i].add (0, seqList.seq (i));
			}

		// Compute initial distance matrix.
		double[][] D = new double [N+1][];
		for (int i = 0; i < N-1; ++ i)
			{
			DnaSequence seq_i = seqList.seq(i);
			double[] D_i = D[i];
			for (int j = i+1; j < N; ++ j)
				{
				DnaSequence seq_j = seqList.seq(j);
				double D_i_j = dcalc.distance (seq_i, seq_j);
				D_i[j] = D_i_j;
				D[j][i] = D_i_j;
				}
			}

		// Set up array of group sizes n_i.
		int[] n = new int [N+1];
		for (int i = 0; i < N; ++ i)
			{
			n[i] = 1;
			}

		//*DEBUG*/ dump ("INITIAL DISTANCE MATRIX", D, n, N);

		// Join trees until only one is left.
		while (N > 1)
			{
			// Find i and j for which D_i_j is smallest.
			double min_d = Double.PositiveInfinity;
			int min_i = 0;
			int min_j = 0;
			for (int i = 0; i < N-1; ++ i)
				{
				double[] D_i = D[i];
				for (int j = i+1; j < N; ++ j)
					{
					double d = D_i[j];
					if (d < min_d)
						{
						min_d = d;
						min_i = i;
						min_j = j;
						}
					}
				}

			//*DEBUG*/ System.out.println ("min_i="+min_i+", min_j="+min_j);

			// Compute node heights for trees <min_i> and <min_j>. Store as
			// branch lengths for now.
			DnaSequenceTree tree_i = tree[min_i];
			tree_i.branchLength (tree_i.root(), 0.5*min_d);
			DnaSequenceTree tree_j = tree[min_j];
			tree_j.branchLength (tree_j.root(), 0.5*min_d);

			// Join trees <min_i> and <min_j>. Add new tree to end of list.
			DnaSequenceTree newtree =
				new DnaSequenceTree (tree_i.length() + tree_j.length() + 1);
			newtree.join (tree_i, tree_j);
			tree[N] = newtree;
			int newn = n[min_i] + n[min_j];
			n[N] = newn;

			// Compute distance from new tree to every other tree.
			double w_i = ((double) n[min_i])/newn;
			double w_j = ((double) n[min_j])/newn;
			double[] D_n = D[N];
			for (int k = 0; k < N; ++ k)
				{
				double[] D_k = D[k];
				double D_n_k = w_i*D_k[min_i] + w_j*D_k[min_j];
				D_n[k] = D_n_k;
				D_k[N] = D_n_k;
				}
			//*DEBUG*/ dump ("DISTANCE MATRIX WITH NEW GROUP ADDED", D, n, N+1);

			// Swap row <N> with row <min_i> and swap row <N-1> with row
			// <min_j>, thus removing rows <min_i> and <min_j> from D, tree, and
			// n.
			double[] swap1 = D[min_i];
			D[min_i] = D[N];
			D[N] = swap1;
			swap1 = D[min_j];
			D[min_j] = D[N-1];
			D[N-1] = swap1;
			DnaSequenceTree swap2 = tree[min_i];
			tree[min_i] = tree[N];
			tree[N] = swap2;
			swap2 = tree[min_j];
			tree[min_j] = tree[N-1];
			tree[N-1] = swap2;
			int swap3 = n[min_i];
			n[min_i] = n[N];
			n[N] = swap3;
			swap3 = n[min_j];
			n[min_j] = n[N-1];
			n[N-1] = swap3;

			// Swap column <N> with column <min_i> and swap column <N-1> with
			// column <min_j>, thus removing columns <min_i> and <min_j> from D.
			for (int i = 0; i <= N; ++ i)
				{
				double[] D_i = D[i];
				double swap4 = D_i[min_i];
				D_i[min_i] = D_i[N];
				D_i[N] = swap4;
				swap4 = D_i[min_j];
				D_i[min_j] = D_i[N-1];
				D_i[N-1] = swap4;
				}

			// Took away two trees and added one.
			D[N] = null;
			tree[N] = null;
			-- N;
			//*DEBUG*/ dump ("DISTANCE MATRIX WITH OLD GROUPS REMOVED", D, n, N);
			}

		// Convert node heights to branch lengths.
		computeBranchLengths (tree[0], tree[0].root());

		// Return the last tree left standing.
		return tree[0];
		}

        /**
	 * Convert the node height stored in the given node of the given tree to a
	 * branch length.
	 *
	 * @param  tree   Tree.
	 * @param  index  Node index.
	 *
	 * @return  Original node height of the given node.
	 */
        private static double computeBranchLengths (DnaSequenceTree tree,int index)
        {
            // Stop recursion when we reach a tip.
            if (index == -1) return 0.0;

            // Convert both children.
            double childHeight = computeBranchLengths(tree, tree.child1(index));
            computeBranchLengths(tree, tree.child2(index));

            // Convert this node.
            Double nodeHeight = tree.branchLength(index);
            if (nodeHeight == null)
            {
                return 0.0;
            }
            else
            {
                tree.branchLength(index, nodeHeight - childHeight);
                return nodeHeight;
            }
        }




    }
}
