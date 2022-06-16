# TreeCont
Scripts for phylogenetic trees and their contractions

Basic scripts for creating a picture of a phylogenetic tree with specified contractions (based on depth and based on node labels). It was tailored for a tree with low thousands of leaves.

# Instalation

Scripts are tailored for Python 3.8, the required libraries can be found in PYTHON_REQUIREMENTS.txt

# Usage

**draw_graphviz_tree.py** generates a dot file from a tree in Stockholm format. You can use DOT (graphviz) to generate a picture. The program also produces HTML files containing a list of subtree labels you can access by clicking on the subtree node.

**draw_tikz_tree.py** generates a tex file from a tree in Stockholm format. You can use LATEX or PDFLATEX to get a PDF. The format here is restricted to A4 paper. The program produces a main tree with contracted nodes, then a subtree for every contracted node. Each node in the main tree is linked to a subtree by a hypertext link. However, only the contracted tree can be printed as well (see preview).

If you want additional information on the parameters, start the script with ```--help``` parameter.

In the directory preview, you can find pdf examples. These were the commands used to produces them:

```
python3 draw_graphviz_tree.py --tree_file tree.sto --dot_name tree_dot_example.dot --do_contraction --helper_labels --dfs_depth 6 --to_contract 3,326,357; dot -Tpdf tree_dot_example.dot >tree_dot_example.pdf
```

HTML files not included in preview.

```
python3 draw_tikz_tree.py --tree_file tree.sto --tex_name tree_tex_example.tex --do_contraction --dfs_depth 18 --to_contract 1,322,358,2794,2791,377,437,447,2774,564,572,1976,2086,2757,2115,2217,2383,2474,1921,1059,1073,1342,582,2925,2872 --only_picture; pdflatex tree_tex_example.tex
```
You can find these files in the **preview** directory.

Feel free to clone and customize.
