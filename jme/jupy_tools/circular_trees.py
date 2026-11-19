"""
Circular (polar) phylogenetic trees with rings of genome metadata, as used
for single-gene and concatenated-marker trees.

Moved from jme.jupy_tools.experimental.vlp_plots_2. The usual entry point is
SyntenyTreePlotter.draw_circular_tree (in tree_plots), which uses the
genome metadata and plot config to build the rings. The pieces here can also
be used directly::

    ring_order = get_ring_order(tree.root, collapsed_dict)
    leaf_arc = 355 / len(ring_order)
    r = draw_circular_tree_recursive(ax, tree.root, leaf_arc=leaf_arc,
                                     collapsed_dict=collapsed_dict)
    for md in tree_metadata_list:   # TreeMetadata objects (plot_metadata)
        r += draw_md_ring(ax, ring_order, md, radius=r, width=.1,
                          leaf_arc=leaf_arc)

Styling:
 * bg_styles_dict / line_styles_dict: {clade: patch/line kwargs}
   get_clade_style_by_type() fills one of these with a background for every
   maximal clade whose leaves are all one type.
 * collapsed_dict: {clade or clade name: True | width | dict(width=..,
   **patch_kwargs)}; get_collapsed_clades_source_host() builds one.
 * support_md: a SupportMetadata to colour/label branch support values.
"""
from collections import defaultdict
from itertools import chain

import numpy
import pandas
from matplotlib import collections, patches
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba

def plot_polar_in_cart(ax, r, th, offset=[0,0], from_degrees=False, **kwargs):
    """ convert r,theta to x,y and plot on cartesian axis """
    x, y = polar_to_cart(r, th, from_degrees)
    ax.plot(x + offset[0],
            y + offset[1],
            **kwargs)

def polar_to_cart(r, th, from_degrees=False):
    if from_degrees:
        th = numpy.radians(th)
        return  (r * numpy.cos(th),
                 r * numpy.sin(th))
    return  (r * numpy.cos(th),
             r * numpy.sin(th))

def get_clade_style_by_type(clade, style_dict, get_leaf_type, type_styles,
                            ignore_nas=False):
    """
    Depth-first crawl of a tree to find clades of all the same type
    PArams:
        node: Bio.Phylo clade
        style_dict: empty dict to save styles to
        get_leaf_type: callable on a leaf node, returns a node type
        type_styles: dict-like mapping types to matplotlib patch kwargs dict
    """
    node_types = set()

    # jest get the type of a leaf node
    if clade.is_terminal():
        return {get_leaf_type(clade.name),}

    # get types of all children
    children_with_one_type = {}
    for child in clade.clades:
        child_types = get_clade_style_by_type(child, style_dict, get_leaf_type,
                                              type_styles, ignore_nas)
        if ignore_nas:
            child_types = set(t for t in child_types if pandas.notna(t))
        if len(child_types) == 1:
            children_with_one_type[child] = next(iter(child_types))
        node_types.update(child_types)

    # if this node has multiple types, set the the style for any children with one type
    if len(node_types) > 1:
        for child, child_type in children_with_one_type.items():
            style_dict[child] = type_styles[child_type]

    return node_types

def get_ring_order(clade, collapsed_dict=None, min_collapsed_size=2):
    if clade in collapsed_dict or clade.name in collapsed_dict:
        N = int(_get_clade_size(collapsed_dict.get(clade,
                                                   collapsed_dict.get(clade.name)),
                               min_collapsed_size)[0])
        return [clade,]*N
    if clade.is_terminal():
        return [clade,]
    return list(chain(*(get_ring_order(child, collapsed_dict,
                                       min_collapsed_size)
                        for child in clade.clades)))


def get_collapsed_clades_source_host(
    tree,
    get_leaf_type,
    known_hosts,
    source_colors,
):
    """
    returns collapsed clades colored by source, but split by host
    """
    collapsed_leaves = set()
    collapsed_clades = {}
    for clade in tree.get_nonterminals():
        consensus = None
        for leaf in clade.get_terminals():
            if leaf in collapsed_leaves:
                # skip if parent clade is collapsed
                break

            source = get_leaf_type(leaf.name)
            host = known_hosts.get(leaf.name, None)

            # only collapse if same source and host (or host is None)
            if consensus is None:
                consensus = (source, host)
            elif consensus[1] is None:
                if source == consensus[0]:
                    consensus = (source, host)
                else:
                    break
            elif host is None:
                if consensus[0] != source:
                    break
            elif consensus != (source, host):
                break
        else:
            # if we get here, ccollapse the clade
            color = source_colors[consensus[0]]
            collapsed_clades[clade] = dict(
                width=clade.count_terminals(),
                ec=color,
                fc=to_rgba(color, .5),
                lw=.5,
            )
            collapsed_leaves.update(clade.get_terminals())
    return collapsed_clades

def _check_triangle(v0, v1, v2):
    """
        # check if the triangle is inverted
        #  IE does the line between v1 and v2 pass between v0 and the origin?
    """
    # the slope of the line from the origin to v0 (its intercept is 0)
    slope_0 = v0[1] / v0[0]
    # The slop of the line between v1 and v2
    slope_12 = (v1[1] - v2[1])/(v1[0] - v2[0])
    # its intercept
    intercept_12 = v1[1] - slope_12 * v1[0]
    # where the 2 lines intercet
    # m1x + b1 = m2x + b2
    # (m1 - m2)x = b2 - b1
    x_intersect = intercept_12 / (slope_0 - slope_12)

    return abs(x_intersect) > abs(v0[0])


def _get_clade_size(collapsed_data, min_collapsed_size):
    """
    Pull the clade size from the collapsed data
    
    collapsed_data is one of:
        int or float indicating size
        dict with a numeric 'width' entry

    any remaining dit entries are return as collapsed_style:
        return clade_size, collapsed_style
    """
    # collapsed value could be one of: True, an integer, or a style
    if isinstance(collapsed_data, bool):
        clade_size = min_collapsed_size
        collapsed_style = {}
    else:
        try:
            clade_size = max(float(collapsed_data), min_collapsed_size)
            collapsed_style = {}
        except (ValueError, TypeError):
            # assume it is a dict
            collapsed_style = dict(collapsed_data)
            # pull width from dict (if it's there)
            try:
                clade_size = max(collapsed_style.pop('width'), min_collapsed_size)
            except KeyError:
                clade_size = min_collapsed_size
    return clade_size, collapsed_style

def _get_clade_sizes(clade, collapsed_dict=None, min_collapsed_size=2):
    if clade in collapsed_dict or clade.name in collapsed_dict:
        N = int(_get_clade_size(collapsed_dict.get(clade,
                                                   collapsed_dict.get(clade.name)),
                               min_collapsed_size)[0])
        return {clade:N}
    if clade.is_terminal():
        return {clade:1}
    clade_sizes = {}
    clade_size = 0
    for child in clade.clades:
        clade_sizes.update(
            _get_clade_sizes(child, collapsed_dict,
                                       min_collapsed_size)
        )
        clade_size += clade_sizes[child]
    clade_sizes[clade] = clade_size
    return clade_sizes

default_line_style = dict(color='k', lw=.5)
def _draw_circular_tree_recursive(ax, clade,
                   distance=0,
                   ax_patches=None,
                   last_distance=None,
                   center=180,
                   total_arc_span=340,
                   min_collapsed_size=2,
                   line_style=None,
                   style_dicts=None,
                   line_styles_dict=None,
                   bg_styles_dict=None,
                   collapsed_dict=None,
                   offset=(0,0),
                   max_distance=None,
                   unit_branch_lengths=None,
                   leaf_arc=None,
                   debug=False,
                   zorder=1,
                   support_md=None,
                  ):
    """
    NOTE: if leaf_arc is supplied, the collapsed_dict muct be keyed on clade
    objects, not names!
    """

    # some things get set automatically based on the first (root) clade
    if max_distance is None:
        max_distance = max(clade.depths().values())
        if max_distance == 0:
            max_distance = max(clade.depths(unit_branch_lengths=True).values())
            unit_branch_lengths = True
        else:
            unit_branch_lengths = False

    if style_dicts is None:
        style_dicts = defaultdict(dict)
        if line_styles_dict is not None:
            style_dicts['lines'] = line_styles_dict
        if bg_styles_dict is not None:
            style_dicts['bg'] = bg_styles_dict
        if collapsed_dict is not None:
            style_dicts['collapsed'] = collapsed_dict

    if 'clade_sizes' not in style_dicts:
        style_dicts['clade_sizes'] = _get_clade_sizes(clade, 
                                                      style_dicts['collapsed'],
                                                      min_collapsed_size
                                                     )
        _n = len(style_dicts['clade_sizes'])
        if debug:
            print(f"N clade sizes: {_n}\n({clade.count_terminals()})")
    if leaf_arc is None:
        # clean up collapsed dict and calculate arc of each thing
        # how many arc segments total
        total_size = 0 # clade.count_terminals()

        # for collapsed nodes, replace leaf count with collapsed size
        skip_nodes = set()
        for cld in clade.get_nonterminals() + clade.get_terminals():
            # is clade collapsed? get the data
            if cld in style_dicts['collapsed']:
                cdata = style_dicts['collapsed'][cld]
                style_dicts['collapsed'][cld.name] = cdata
            elif cld.name in style_dicts['collapsed']:
                cdata = style_dicts['collapsed'][cld.name]
                # add duplicate entry with clade object
                style_dicts['collapsed'][cld] = cdata
            else:
                # clade is not collapsed
                if cld.is_terminal():
                    total_size += 1
                continue
            # don't consider clades inside other collapsed clades
            if not any(t in skip_nodes for t in cld.get_terminals()):
                c_size, _ = _get_clade_size(cdata, min_collapsed_size)
                total_size += c_size
                skip_nodes.update(cld.get_terminals())

        # get fraction of requested span
        leaf_arc = total_arc_span / total_size

    if unit_branch_lengths:
        clade_depths = clade.depths(unit_branch_lengths=True)
        def get_clade_distance(node):
            return clade_depths[node]
    else:
        get_clade_distance = clade.distance

    if ax_patches is None:
        add_patches_on_exit = True
        ax_patches = defaultdict(list)
    else:
        add_patches_on_exit = False

    # get the line style. Use clade specific if available, fall back to parent style, then to default
    if clade in style_dicts['lines']:
        # keep default settings for anything not overriden by the clade style
        line_style = dict(default_line_style)
        line_style.update(style_dicts['lines'][clade])
    elif line_style is None:
        line_style = default_line_style

    # draw radial line from parent arc
    if last_distance is not None:
        tot_distance = last_distance + distance
        plot_polar_in_cart(ax,
                           (last_distance, tot_distance),
                           (center, center),
                           offset=offset,
                           from_degrees=True,
                           zorder=zorder,
                           **line_style)
    else:
        tot_distance = distance

    # is this clade collapsed?
    collapsed = style_dicts['collapsed'].get(clade, None)

    # figure out where this clade is
    N = style_dicts['clade_sizes'][clade]
    start = center - leaf_arc * (N/2)
    end = center + leaf_arc * (N/2)
    if debug:
        print((clade.count_terminals() if clade.name is None else clade.name), N, start, center, end)

    # if clade has a background set
    bg_style = style_dicts['bg'].get(clade, None)
    if bg_style is not None:
        ax_patches['wedges'].append(patches.Wedge(offset,
                                   r=tot_distance,
                                   theta1=start,
                                   theta2=end,
                                   width=tot_distance - max_distance,
                                   **bg_style,
                                  )
                    )

    # Yes if there is a nonnull or False value. Empty dict is OK (but would resolve to False)
    if collapsed is not None and collapsed is not False:
        if debug:
            print(f"collapsing node ({collapsed})")

        triangle_width, collapsed_style_data = \
            _get_clade_size(collapsed, min_collapsed_size)

        # otherwise, use the clade size as the width
        # get default style
        # use bg settings with line style overrides
        triangle_style = dict(line_style)
        if bg_style is not None:
            triangle_style.update(bg_style)

        # we pulled any style data from the collapsed dict earlier.
        # Add it here:
        triangle_style.update(collapsed_style_data)

        # get triangle shape
        distances = numpy.array([get_clade_distance(t) for t in clade.get_terminals()])
        if debug:
            print(f"{len(distances)} distances from {distances.min()} to "
                  f"{distances.max()}")
        # set up triangle in polar corrdinates
        triangle_verts_polar_unzipped = [
            (tot_distance, tot_distance + distances.max(), tot_distance + distances.min()),
            (center, start, end)
        ]
        if debug:
            print(f"Trinagle locs: {triangle_verts_polar_unzipped}")
        # convert to cartesian pairs
        triangle_verts = list(zip(*polar_to_cart(*triangle_verts_polar_unzipped,
                                                 from_degrees=True)))
        if debug:
            print(f"Converted trinagle locs: {triangle_verts}")

        if _check_triangle(*triangle_verts):
            ax.add_patch(patches.Polygon(triangle_verts,
                                         zorder=zorder,
                                         **triangle_style),
                        )
        else:
            ax.add_patch(
                    patches.Arc(
                        offset,
                        width=tot_distance*2,
                        height=tot_distance*2,
                        theta1=start,
                        theta2=end,
                        zorder=zorder,
                        **line_style
                    )
            )
        depths = distances + tot_distance
    else:
        # figure out where children go and draw them
        last_end = start
        centers = []
        depths = [tot_distance,]
        for child in clade.clades:
            d = get_clade_distance(child)    # radial distance
            clade_size = style_dicts['clade_sizes'][child]
            w2 = clade_size * leaf_arc / 2        # 1/2 angle covered
            c = last_end + w2                     # center (half from previous end)
            centers.append(c)
            d = _draw_circular_tree_recursive(
                               ax, child,
                               ax_patches=ax_patches,
                               distance=d,
                               last_distance=tot_distance,
                               center=c,
                               offset=offset,
                               line_style=line_style,
                               style_dicts=style_dicts,
                               leaf_arc=leaf_arc,
                               max_distance=max_distance,
                               unit_branch_lengths=unit_branch_lengths,
                               zorder=zorder,
                               support_md=support_md,
                               debug=debug,
                              )
            depths.append(d)
            last_end = c + w2

        # connect children with arc
        if centers:
            # if it's not the first branch
            if tot_distance > 0:
                diameter = tot_distance * 2
                arc = \
                    patches.Arc(
                        offset,
                        width=diameter, height=diameter,
                        theta1=min(centers),
                        theta2=max(centers),
                        zorder=zorder,
                        **line_style
                    )
                #ax_patches['arcs'].append(arc)
                ax.add_patch(arc)

            # if there are children and a confidence value
            if clade.confidence is not None:
                if support_md:
                    # color the branch point by confidence
                    annotate_clade_support(ax,
                                           clade.confidence, 
                                           center,
                                           tot_distance,
                                           support_md)


    if add_patches_on_exit:
        ax.add_collection(collections.PatchCollection(ax_patches['wedges'],
                                                      zorder=zorder,
                                                      match_original=True),
                         )

    return max(depths)


def annotate_clade_support(ax, confidence, theta, radius, support_md):
    """
    draw a circle at (theta, radius) colored by confidence and support_md
    """
    x, y = polar_to_cart(radius, theta, from_degrees=True)
    ax.scatter(x, y, color=support_md.get_color(confidence), **support_md.plot_args)

class SupportMetadata():
    def __init__(self, **kwargs):
        self.color_params = kwargs.pop('color_params', {})
        self.color_map = plt.get_cmap(self.color_params.get('color_map',
                                                            'binary'))
        self.plot_args = kwargs

    def get_color(self, confidence):
        return self.color_map(confidence)


# Set up collapsed clades

def draw_md_ring(ax, ring_order, mdmd, radius=3, width=.1, center=180, leaf_arc=None, offset=[0,0]):
    width = width * mdmd.relative_size
    N = len(ring_order)
    if leaf_arc is None:
        leaf_arc = 340 / N

    ax_patches = []
    start = center - leaf_arc * (N/2)
    for leaf in ring_order:
        color = mdmd.get_item_color(leaf)
        if mdmd.size_norm is not None:
            leaf_size = mdmd.get_item_size(leaf) * width
        else:
            leaf_size = width
        end = start + leaf_arc

        if True:
            name = (leaf.name if leaf.is_terminal() else
            f"{leaf.count_terminals()} leaves")

        ax_patches.append(patches.Wedge(offset,
                                   r=radius,
                                   theta1=start,
                                   theta2=end,
                                   width=-leaf_size,
                                   fc=color,
                                   ec=None,
                                  )
            )
        start = end

    ax.add_collection(collections.PatchCollection(ax_patches, match_original=True))
    return width


# public name for the recursive drawer
draw_circular_tree_recursive = _draw_circular_tree_recursive
