
process hashing_concat{
   
    cache 'lenient'
    debug true

    input:
    path hashing_demux_anndata

    output:
    path "concatenated_hashing_demux.h5ad", emit: concatenated_hashing_demux

    script:
    """
    hashing_concat.py ${hashing_demux_anndata} --output concatenated_hashing_demux.h5ad
    """
}
