"""
Functions for drawing gee-centric genome maps



"""

def draw_genes(gene_data,
               location=0,
               thickness=1,
               transform=None,
               orientation='horizontal',
              ):
    """
    params:
        gene_data:
            DataFrame of annotations with columns:
                start
                end
                strand (optional)
                label (optional)
                color (optional)
        location:
            y coordinate to draw on (or x if orientation is 'vertical')
        thickness:
            vertical thickness of arrows (or horizontal if orientation
                                          is set to 'vertical')
        transform:
            transformation to apply to coordinates
        orientation:
            defaults to horizontal
    """

