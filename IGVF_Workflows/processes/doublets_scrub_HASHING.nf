
process doublets_scrub_HASHING {
 
    cache 'lenient'
    
    input:
        path mudata

    output:
        path "mdata_doublets.h5mu", emit: mudata_doublet

    script:
        """
        doublets_HASHING.py --input ${mudata}
        """
}
