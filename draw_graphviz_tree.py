import argparse
import math
import networkx as nx
from Bio import Phylo as ph
import os

parser = argparse.ArgumentParser()
parser.add_argument("--tree_file",
                    default="tree.fst", help="Path to the input tree in Stockholm format.",
                    type=str)
parser.add_argument("--entry_width",
                    default=15, help="Width od the single graphviz node",
                    type=int)
parser.add_argument("--dot_name",
                    default="testing_dot.dot", help="Path to the output dot file.",
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
parser.add_argument("--helper_labels",
                    default=False,
                    dest='helper_labels', help="Print unique node identifiers. "
                                               "Handy for picking specific trees to contract.",
                    action='store_true')
parser.add_argument("--entry_height",
                    default=10, help="Heigth od the single graphviz node",
                    type=int)
parser.add_argument("--fontsize",
                    default=50, help="Fontsize of graphviz node label.",
                    type=int)


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


def main(args):
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

    def build_subtree_label(node, node_info):
        subtree_size = len([x for x in node_info["name"] if x is not None])  # no of sequences in subtree

        if args.helper_labels:
            label = f"v={node}\n{subtree_size} seqs"
        else:
            label = f"{subtree_size} seqs"

        return label

    def make_html_file(node, node_info):
        dirname = args.dot_name[:-4]
        filename = f"{dirname}/subtree_{node}.html"
        os.makedirs(dirname, exist_ok=True)

        with open(filename, mode='w') as writer:
            print("<!DOCTYPE HTML>\n<html>", file=writer)
            print("<body>", file=writer)
            print("<table>", file=writer)
            # header
            print("<tr>", file=writer)
            print("<th>Sequence identificator</th>\n", file=writer)
            print("</tr>", file=writer)
            # elements

            for name in node_info["name"]:
                if name is None:
                    continue

                colors = node_info["color"]
                if len(colors) == 0:
                    print(f"<tr style=\"background-color: {select_color(colors)}\">", file=writer)
                else:
                    print(f"<tr>", file=writer)
                print(f"<td>{name}</td>\n",
                      file=writer)
                print("</tr>", file=writer)
            print("</table>", file=writer)
            print("<br>", file=writer)

            print("</body>", file=writer)
            print("</html>", file=writer)

        return filename

    def select_color(color_options):
        return color_options.pop()

    def build_leaf_label(node, node_info, row=30):
        name = node_info["name"]
        t = '\n'.join([name[x * row:x * row + row] for x in range(math.ceil(len(name) / row))])
        return t

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

    new_graph, info_on_node = load_tree(args.tree_file)

    if args.do_contraction:
        if args.to_contract is not None:
            node_ids_to_contract = args.to_contract.split(',')
        else:
            node_ids_to_contract = args.to_contract
        do_contraction(node_ids_to_contract, args.dfs_depth)

    with open(args.dot_name, mode='w') as treewr:
        print("graph {", file=treewr)
        # print("rankdir=\"LR\"", file=treewr)
        print("splines=\"false\"", file=treewr)
        print("overlap=\"false\"", file=treewr)
        print(f"ranksep=\"{args.entry_width // 2}\"", file=treewr)
        print("nodesep=\"2\"", file=treewr)
        print(f"fontsize=\"{args.fontsize}\"", file=treewr)

        for vertex in new_graph.nodes:
            node_info = info_on_node[vertex]
            if node_info["contracted_sequences"] is not None:
                label = build_subtree_label(vertex, node_info)
                htmlfile = make_html_file(vertex, node_info)

                if len(node_info["color"]) == 0:
                    color = args.neutral_color
                elif len(node_info["color"]) == 1:
                    color = node_info["color"].pop()
                else:
                    color = select_color(node_info["color"])
                print(
                    f"{vertex} [shape=\"triangle\", color=\"black\", width={args.entry_width}, style=\"filled\", "
                    f"fillcolor=\"{color}\", label=\"{label}\", fontsize={args.fontsize}, URL=\"{htmlfile}\"]",
                    file=treewr)
            else:
                represented_clade = node_info["clade_info"]
                if represented_clade.name is None:
                    if represented_clade.confidence is None:
                        print(f"{vertex} [shape=\"point\", color=\"black\"]",
                              file=treewr)
                    else:
                        conf = "{:.2f}".format(represented_clade.confidence)
                        if args.helper_labels:
                            print(
                                f"{vertex} [shape=\"box\", width=1, color=\"black\", label=\"v={vertex}\", "
                                f"fontsize={args.fontsize}]",
                                file=treewr)
                        else:
                            print(f"{vertex} [shape=\"point\", color=\"black\"]",
                                  file=treewr)
                else:
                    label = build_leaf_label(vertex, node_info)
                    color = select_color(node_info["color"])
                    print(
                        f"{vertex} [shape=\"box\", color=\"black\", width={args.entry_width}, style=\"filled\", "
                        f"fillcolor=\"{color}\", label=\"{label}\", fontsize={args.fontsize}]",
                        file=treewr)

        for u, v in new_graph.edges:
            # print(f"{u} -- {v} [headport=w, tailport=e];", file=treewr)
            print(f"{u} -- {v} [headport=n, tailport=s];", file=treewr)

        print("}", file=treewr)


if __name__ == '__main__':
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
