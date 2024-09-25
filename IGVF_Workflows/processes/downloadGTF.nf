
process downloadGTF {
    input:
    val gtf_url

    output:
    path "gencode_gtf.gtf.gz" , emit: gencode_gtf

    script:
    """
        wget "${gtf_url}" -O gencode_gtf.gtf.gz
    """
}
