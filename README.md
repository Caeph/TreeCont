# TreeCont
Scripts for phylogenetic trees and their contractions

Basic scripts for creating a picture of a phylogenetic tree with specified contractions (based on depth and based on node labels). It was tailored for a tree with low thousands of leaves.

# Instalation

Scripts are tailored for Python 3.8, the required libraries can be found in PYTHON_REQUIREMENTS.txt

# Usage

**draw_graphviz_tree.py** generates a dot file from a tree in Stockholm format. You can use DOT (graphviz) to generate a picture. The program also produces HTML files containing a list of subtree labels you can access by clicking on the subtree node.

**draw_tikz_tree.py** generates a tex file from a tree in Stockholm format. You can use LATEX or PDFLATEX to get a PDF. The format here is restricted to A4 paper. The program produces a main tree with contracted nodes, then a subtree for every contracted node. Each node in the main tree is linked to a subtree by a hypertext link.

In the directory preview, you can find pdf examples. These were the commands used to produces them:

Feel free to clone and customize.
