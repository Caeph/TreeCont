import networkx as nx
from Bio import Phylo as ph
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tree_file",
                    default="tree.fst",
                    type=str,
                    help="Path to the input tree in Stockholm format."
                    )
parser.add_argument("--tex_name",
                    default="testing_dot.tex", help="Path to the output tex file.",
                    type=str)
parser.add_argument("--to_contract",
                    default=None, help="Labels of nodes to contract divided by commas (no spaces). You can get "
                                       "the labels using the graphviz drawing "
                                       "software with --helper_labels. Labels are interchangable. Example:"
                                       "1,2,3,58,615",
                    type=str)
parser.add_argument("--do_contraction",
                    default=False,
                    dest='do_contraction', help="Contraction switch. If not specified, contraction is not performed",
                    action='store_true')
parser.add_argument("--dfs_depth",
                    default=5,
                    type=int,
                    help="Maximal depth of the tree. If not specified, dfs contraction is not performed. Must be used "
                         "with the <do_contraction> switch."
                    )
parser.add_argument("--only_picture",
                    default=False,
                    dest='only_picture', action='store_true', help="Print only the big contracted picture, "
                                                                   "do not draw the subtrees."
                    )

basic_conf_col = "black"
confidence_colors = [
    (0.9, f"{basic_conf_col}"),
    (0.7, f"{basic_conf_col}!50"),
    (0.5, f"{basic_conf_col}!10"),
]
low_conf = 0.5

title = "A SUPER COOL PHYLOGENETIC TREE"


def print_preabmle(treewr, TITLE, args):
    print("\\documentclass{article}", file=treewr)
    print("\\usepackage[x11names, svgnames, rgb]{xcolor}", file=treewr)
    print("\\usepackage[utf8]{inputenc}", file=treewr)
    print("\\usepackage{tikz}", file=treewr)
    print("\\usepackage{scalefnt}", file=treewr)
    print("\\usepackage{longtable}", file=treewr)
    print("\\usetikzlibrary{snakes,arrows,shapes,calc}", file=treewr)
    print("\\usepackage{amsmath}", file=treewr)
    print("\\usepackage{colortbl}", file=treewr)
    print("\\usepackage{geometry}", file=treewr)
    print("\\geometry{a4paper, total={170mm,257mm}, left=10mm, top=5mm}", file=treewr)
    print("\\usepackage[hidelinks]{hyperref}", file=treewr)

    # use arial font
    print("\\usepackage{helvet}\n\\renewcommand{\\familydefault}{\\sfdefault}", file=treewr)

    print("\\begin{document}\n\\pagestyle{empty}\n\\enlargethispage{100cm}", file=treewr)

    print("\\tikzset{\nhyperlink node/.style={\n"
          "alias=sourcenode,\n"
          "append after command={\n"
          "let     \\p1 = (sourcenode.north west),\n"
          "\\p2=(sourcenode.south east),"
          "\\n1={\\x2-\\x1},\n"
          "\\n2={\\y1-\\y2} in\n"
          "node [inner sep=0pt, outer sep=0pt,anchor=north west,at=(\\p1)] "
          "{\\hyperlink{#1}{\\XeTeXLinkBox{\\phantom{\\rule{\\n1}{\\n2}}}}}\n"
          "}\n}\n}", file=treewr)
    print("\n", file=treewr)

    if not args.only_picture:
        print("\\vspace*{8cm}", file=treewr)
        print("{\\centering\\Large\\bfseries " + TITLE + "}", file=treewr)
        print("\\\\", file=treewr)
        # print("\\vfill{}", file=treewr)
        print("{\\centering Click the taxon to get to the subtree details.}", file=treewr)
        print("\\newpage", file=treewr)


def generate_tree(treewr, treegraph, root='0', LR=False,
                  paperwidth=257,
                  nodestyle=lambda node: "draw=black,fill,rectangle,minimum height=10mm,minimum width=0.01cm",
                  labelgen=lambda node: "",
                  anchor="south east",
                  special_features=lambda node, x, y, order, maxorder: None,
                  paperheight=190, heightstep=None, widthstep=None
                  ):
    # kazdemu nodu urcit heightlvl a poradi v ranku
    ranks_height = {}
    ranks_width = {}

    preds = nx.dfs_predecessors(treegraph, root)
    succs = nx.dfs_successors(treegraph, root)

    for vertex in treegraph.nodes():
        current = vertex
        level = 0
        while current != root:
            level += 1
            current = preds[current]
        ranks_height[vertex] = level

    postorder = nx.dfs_postorder_nodes(treegraph, source=root)
    leaves = [n for n, deg in treegraph.out_degree() if deg == 0]

    widthlevel = 0
    for node in postorder:
        if node in leaves:
            ranks_width[node] = widthlevel
            widthlevel += 1

    notleafs = sorted([n for n, deg in treegraph.out_degree() if deg > 0],
                      key=lambda n: -ranks_height[n])
    for notleaf in notleafs:
        succ_widths = [ranks_width[s] for s in succs[notleaf]]
        ranks_width[notleaf] = sum(succ_widths) / len(succ_widths)

    # from nodes counts in levels assign step sizes if undefined
    if widthstep is None:
        widthstep = paperwidth / max(ranks_width.values())
    if heightstep is None:
        heightstep = paperheight / max(ranks_height.values())

    print("\\begin{center}\n\\begin{tikzpicture}[>=latex',line join=bevel,"
          "cross/.style={path picture={ \\draw[black] (path picture bounding box.south east) -- "
          "(path picture bounding box.north west) (path picture bounding box.south west) -- "
          "(path picture bounding box.north east);}}]"
          "]", file=treewr)

    if not LR:
        maxorder = max(ranks_height.values())
    else:
        maxorder = max(ranks_width.values())

    # print vertices
    positions = {}  # vertex : "xmm,ymm"
    for vertex in treegraph.nodes():
        height, width = ranks_height[vertex], ranks_width[vertex]
        nodelook = nodestyle(vertex)
        if not LR:
            x = width * widthstep
            y = height * heightstep
            order = height
        else:
            x = height * heightstep
            y = width * widthstep
            order = width
        positioning = f"{x}mm,{y}mm"
        positions[vertex] = (x, y)
        print(f"\\node ({vertex}) at ({positioning}) [{nodelook}]" + " {};", file=treewr)
        # node font=\\tiny

        additional = special_features(vertex, x, y, order, maxorder)
        if additional is not None:
            print(additional, file=treewr)

    # rooting
    rootx, rooty = positions[root]
    print(f"\\node (root) at ({rootx - heightstep}mm,{rooty}mm) [draw=black,ultra thin,fill,text width=0.01mm,"
          f"inner sep=0pt,rectangle]" + " {};", file=treewr)
    print(f"\\draw [very thick] (root) -- ({root});", file=treewr)

    # print edges
    for from_, to_ in treegraph.edges():
        posi_from = positions[from_]
        posi_to = positions[to_]
        if LR:
            midpoint = f"{posi_from[0]}mm,{posi_to[1]}mm"
        else:
            midpoint = f"{posi_to[0]}mm,{posi_from[1]}mm"
        print(f"\\node ({from_}{to_}midpoint) at ({midpoint}) [draw=black,ultra thin,fill,text width=0.03mm,"
              f"inner sep=0pt,rectangle]" + " {};", file=treewr)
        print(f"\\draw [very thick] ({from_}) -- ({midpoint});", file=treewr)
        print(f"\\draw [very thick] ({midpoint}) -- ({to_});", file=treewr)

    for vertex in treegraph.nodes():
        print(f"\\node ({vertex}label) at ({vertex}.west) [anchor={anchor}, font=\\tiny]" + " {" + labelgen(
            vertex) + "};", file=treewr)

    print("\\end{tikzpicture}\n\\end{center}", file=treewr)


def main(args):
    def load_tree(tree_file):
        tree = ph.read(tree_file, 'newick')
        graph = ph.to_networkx(tree)
        node_info = {}

        clade_nodes = {}

        new_graph = nx.DiGraph()

        for i, node in enumerate(graph.nodes):
            clade_nodes[node] = i
            name = node.name
            node_info[f"{i}"] = {
                "clade_info": node,
                "name": name,
                "contracted_sequences": None,
                "color": ["grey"]  # can be adjusted through an external dictionary
            }

        for a, b in graph.edges:
            current = ["", ""]
            for i, node in enumerate([a, b]):
                current[i] = f"{clade_nodes[node]}"

            for item in current:
                new_graph.add_node(item)
            new_graph.add_edge(*current)

        return new_graph, node_info

    def do_contraction(to_contract, dfs_depth):
        # contract specified nodes
        if to_contract is not None:
            for contract_node in to_contract:
                node_contraction(contract_node)

        if dfs_depth is not None:
            shallow_tree = nx.dfs_tree(new_graph, '0', dfs_depth)
            dfs_to_contract = ([n for n, deg in shallow_tree.out_degree() if deg == 0])
            for contract_node in dfs_to_contract:
                node_contraction(contract_node)

    def node_contraction(contract_node):
        if new_graph.out_degree(contract_node) == 0:
            return

        subtree = nx.dfs_tree(new_graph, contract_node)
        following_nodes = subtree.nodes()

        # contract info on the node
        clade_infos, names, color = [], [], []
        for x in following_nodes:
            if info_on_node[x]["contracted_sequences"]:  # is contracted
                clade_infos.extend(info_on_node[x]["clade_info"])
                names.extend(info_on_node[x]["name"])

                # color choosing -- this is fairly stupid and unusable
                color.extend(info_on_node[x]["color"])

            else:  # is not contracted
                clade_infos.append(info_on_node[x]["clade_info"])
                names.append(info_on_node[x]["name"])
                color = info_on_node[x]["color"]

        # manipulate_graph, save contracted info
        new_graph.remove_nodes_from([x for x in following_nodes if x != contract_node])
        info_on_node[contract_node] = {
            "clade_info": clade_infos,
            "contracted_sequences": True,
            "name": names,
            "color": color
        }

    new_graph, info_on_node = load_tree(args.tree_file)
    entire_tree = new_graph.copy()
    _, entire_tree_info_on_node = load_tree(args.tree_file)

    if args.do_contraction:
        node_ids_to_contract = [x for x in args.to_contract.split(',') if len(x) > 0]
        do_contraction(node_ids_to_contract, args.dfs_depth)

    with open(args.tex_name, mode='w') as treewr:
        # preambule
        print_preabmle(treewr, title, args)
        root = '0'

        def style_generator(vertex):
            if new_graph.out_degree(vertex) > 0:
                node_info = info_on_node[vertex]
                represented_clade = node_info["clade_info"]
                if (represented_clade.confidence is None) or (represented_clade.confidence < low_conf):
                    return f"draw=black,ultra thin,fill,rectangle,text width=0.1mm,inner sep=0pt"
                for thr, col in confidence_colors:
                    if represented_clade.confidence >= thr:
                        return f"draw=black, ultra thin, circle, fill={col}, minimum width=2mm, inner sep=0pt"

            return "fill=black, rectangle, minimum height=0.55cm, text width=0.5mm, font={\\tiny}, inner sep=0pt, " \
                   f"hyperlink node=subtree{vertex}"

        def label_generator(vertex):
            if new_graph.out_degree(vertex) > 0:
                return ""
            node_info = info_on_node[vertex]
            if node_info["contracted_sequences"] is None:
                return "1"
            subtree_size = len([x for x in node_info["name"] if x is not None])
            return f"{subtree_size}"

        def style_generator_small(vertex):
            if new_graph.out_degree(vertex) > 0:
                node_info = info_on_node[vertex]
                represented_clade = node_info["clade_info"]
                if (represented_clade.confidence is None) or (represented_clade.confidence < low_conf):
                    return f"draw=black,ultra thin,fill,rectangle,text width=0.1mm,inner sep=0pt"
                for thr, col in confidence_colors:
                    if represented_clade.confidence >= thr:
                        return f"draw=black, ultra thin, circle, fill={col}, minimum width=2mm, inner sep=0pt"

            return f"draw=black, fill, rectangle, minimum height=1cm, text width=0.5mm, inner sep=0pt"

        def style_generator_subtree(vertex):
            node_info = entire_tree_info_on_node[vertex]
            if entire_tree.out_degree(vertex) > 0:
                represented_clade = node_info["clade_info"]
                if (represented_clade.confidence is None) or (represented_clade.confidence < low_conf):
                    return f"draw=black,ultra thin,fill,rectangle,text width=0.5mm,inner sep=0pt"
                for thr, col in confidence_colors:
                    if represented_clade.confidence >= thr:
                        return f"draw=black, ultra thin, circle, fill={col}, minimum width=1mm,inner sep=0pt"
            else:
                return f"draw=black,ultra thin,rectangle,fill=black,anchor=west,text width=1mm,inner sep=0pt"

        # Colapsed overall tree -- generate tikz file
        generate_tree(treewr, new_graph, LR=True,
                      nodestyle=style_generator,
                      labelgen=label_generator,
                      heightstep=7,
                      paperwidth=180)  # 270

        # Subtrees
        if not args.only_picture:
            contracted = [n for n, deg in new_graph.out_degree() if deg == 0]
            for contr_node in contracted:
                print("\\newpage", file=treewr)
                print("\\hypertarget{subtree" + contr_node + "}{\\section{Subtree details}}", file=treewr)
                print("Click the node in the detailed tree to get to the details of the molecule.", file=treewr)
                print("\\hspace*{-0.8cm}", file=treewr)

                # node img in the contracted graph
                shallow_tree = nx.dfs_tree(new_graph, contr_node)
                generate_tree(treewr, shallow_tree, root=contr_node,
                              LR=True,
                              nodestyle=style_generator_small,
                              labelgen=label_generator,
                              widthstep=30,
                              heightstep=7.5, paperwidth=50)

                # entire subtree picture
                subtree = nx.dfs_tree(entire_tree, contr_node)
                w = 1.1

                def moved_labels_subtree(vertex, x, y, order, maxorder):
                    if entire_tree.out_degree(vertex) > 0:
                        return ""
                    node_info = info_on_node[vertex]
                    name = node_info["name"].replace("_", " ")  # cant use _ in texfile

                    name = f"\\tiny {name}"

                    if (order % 2) == 0:
                        result = f"\\node ({vertex}movedlabel) [anchor=west,minimum size=4mm,inner sep=0pt, " \
                                 f"] at ({x + 55}mm,{y}mm) " + "{" + name + "};\n"
                        result += f"\\draw [thin] ({vertex}) -- ({vertex}movedlabel);"
                    else:
                        result = f"\\node({vertex}movedlabel) [anchor=west,minimum size=4mm,inner sep=0pt] at ({x + 3}mm,{y}mm) " + "{" + name + "};"

                    return result

                generate_tree(treewr, subtree, root=contr_node,
                              LR=True,
                              nodestyle=style_generator_subtree,
                              anchor="west",
                              special_features=moved_labels_subtree,
                              widthstep=w,
                              heightstep=3)

                # here, an exhaustive table can be defined
        print("\\end{document}", file=treewr)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
