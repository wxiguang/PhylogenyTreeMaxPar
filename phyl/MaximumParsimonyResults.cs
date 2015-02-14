using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phyl
{
    class MaximumParsimonyResults 
    {

        
        // Extra padding in tree list.
        private static int PAD = 32;
        // Hidden data members.
        // List of tree signatures, its size, and its capacity.
        private int[][] treeList;
        private int size;
        private int capacity;
        // Parsimony score.
        public int score;
        // For detecting modifications during an iteration.
        private int modCount;
        // Extra padding to avert cache interference.
        private static long p0, p1, p2, p3, p4, p5, p6, p7;
        private static long p8, p9, pa, pb, pc, pd, pe, pf;

        /**
        * Construct a new uninitialized maximum parsimony results object. This
        * constructor is for use only by object deserialization.
        */
        public MaximumParsimonyResults()
        {
        }
        /**
        * Construct a new maximum parsimony results object. The tree list is
        * initialized to an empty list with the given capacity. The parsimony score
        * is initialized to <TT>Integer.MAX_VALUE</TT>.
        *
        * @param capacity Capacity.
        *
        * @exception IllegalArgumentException
        * (unchecked exception) Thrown if <TT>capacity</TT> &le; 0.
        */
        public MaximumParsimonyResults (int capacity)
        {
            if (capacity <= 0)
            {
                throw new ArgumentException ("MaximumParsimonyResults(): capacity (= " + capacity + ") illegal");
            }
            this.size = 0;
            this.capacity = capacity;
            this.treeList = new int[capacity + PAD][];
            this.score = int.MaxValue;
        }

        /**
        * Construct a new maximum parsimony results object that is a copy of the
        * given maximum parsimony results object.
        *
        * @param results Maximum parsimony results object to copy.
        */
        public MaximumParsimonyResults (MaximumParsimonyResults results)
        {
            this.size = results.size;
            this.capacity = results.capacity;
            this.treeList = new int[results.capacity + PAD][];
            for (int i = 0; i < results.size; ++i)
            {
                this.treeList[i] = (int[])results.treeList[i].Clone();
            }
            this.score = results.score;
        }





        /**
        * Clear this maximum parsimony results object. Afterwards, the tree list is
        * empty and the parsimony score is <TT>Integer.MAX_VALUE</TT>.
        */
        public void clear()
        {
            ++modCount;
            for (int i = 0; i < size; ++i) treeList[i] = null;
            size = 0;
            score = int.MaxValue;
        }

        /**
        * Add the given tree with the given parsimony score to this maximum
        * parsimony results object. The following invariant is maintained: This
        * maximum parsimony results object contains only those trees with the
        * smallest parsimony score seen so far; and only the first <I>C</I> such
        * trees are stored, where <I>C</I> is the capacity.
        *
        * @param tree Tree signature.
        * @param score Tree's parsimony score.
        */
        public void add (int[] tree,int score)
        {
            ++modCount;
            if (score < this.score)
            {
                clear();
                this.score = score;
            }
            if (score == this.score && size < capacity)
            {
                treeList[size] = (int[])tree.Clone();
                ++size;
            }
        }

        public void addAll(MaximumParsimonyResults results)
        {
            ++modCount;
            if (results.score < this.score)
            {
                clear();
                this.score = results.score;
            }
            if (results.score == this.score)
            {
                int i = 0;
                while (this.size < this.capacity && i < results.size)
                {
                    this.treeList[this.size] = (int[])results.treeList[i].Clone();
                    ++this.size;
                    ++i;
                }
            }
        }

        /**
        * Returns the tree at the given index in this maximum parsimony results
        * object.
        *
        * @param i Index, 0 &le; <TT>i</TT> &le; <TT>size()-1</TT>.
        *
        * @return Tree signature at index <TT>i</TT>.
        *
        * @exception IndexOutOfBoundsException
        * (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
        */
        public int[] tree
        (int i)
        {
            if (0 > i || i >= size)
            {
                throw new IndexOutOfRangeException ("MaximumParsimonyScore.tree(): i (= " + i + ") out of bounds");
            }
            return treeList[i];
        }

        /**
        * Reduce this maximum parsimony results object's score to the given score.
        * If this object's score is less than or equal to <TT>score</TT>, this
        * object is unchanged. If this object's score is greater than
        * <TT>score</TT>, this object is cleared and its score is set to
        * <TT>score</TT>.
        *
        * @param score Parsimony score.
        */
        public void reduceScore
        (int score)
        {
            if (this.score > score)
            {
                clear();
                this.score = score;
            }
        }


    }
}
