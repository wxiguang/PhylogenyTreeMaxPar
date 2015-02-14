using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace phyl
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }


        static int N=5;
        static DnaSequenceTree[] treestack;
        static DnaSequence[] seqstack;


        /**
         * Generate all trees with one more tip node than the tree at the given
         * stack level.
         *
         * @param  level  Stack level.
         */
        private static void generateTrees (int level)
		{
		DnaSequenceTree tree = treestack[level];
		if (level == N-1)
			{


			MessageBox.Show ("Tree No: " + tree.toString());
			}
		else
			{
			DnaSequenceTree nexttree = treestack[level+1];
			int len = tree.length();
			for (int i = 0; i < len; ++ i)
				{
				nexttree.copy (tree);
				nexttree.add (i, seqstack[level+1]);
				generateTrees (level+1);
				}
			}
		}


        	/**
	 * Returns the DNA sequence name determined by the given level. Level 0 =
	 * <TT>"A"</TT>, level 1 = <TT>"B"</TT>, and so on.
	 *
	 * @param  level  Level.
	 *
	 * @return  DNA sequence name.
	 */
	private static String nameForLevel
		(int level)
		{
		StringBuilder buf = new StringBuilder();
		int v = level + 1;
		while (v > 0)
			{
			buf.Insert (0, v2c[v%27]);
			v /= 27;
			}
		return buf.ToString();
		}

	private static char[] v2c = new char[]
		{' ', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
		 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
		 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

        public void inital() {

            DnaSequenceList list = DnaSequenceList.read();
            MessageBox.Show(list.mySequence[3].mySites[11].ToString());

            	int sites = list.seq(0).length();
                MessageBox.Show("Excising uninformative sites ...");
                int nChanges = list.exciseUninformativeSites();
            	MessageBox.Show (sites + " sites");
                MessageBox.Show(list.seq(0).length() + " informative sites");
                MessageBox.Show(nChanges + " state changes from uninformative sites");



        }

        private void button1_Click(object sender, EventArgs e)
        {
            //inital();

            treestack = new DnaSequenceTree[N];
            for (int i = 0; i < N; ++i)
            {
                treestack[i] = new DnaSequenceTree(2 * i + 1);
            }
            seqstack = new DnaSequence[N];
            for (int i = 0; i < N; ++i)
            {
                seqstack[i] = new DnaSequence(0, 0, nameForLevel(i));
            }
            treestack[0].add(0, seqstack[0]);
            generateTrees(0);
        }
    }
}
