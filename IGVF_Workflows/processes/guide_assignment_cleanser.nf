process guide_assignment_cleanser {

    input:
        path mudata_input
        val threshold

    output:
        path "cleanser_mudata_output.h5mu", emit: guide_assignment_mudata_output

    script:
        """
            cleanser.py  -i ${mudata_input} -o cleanser_mudata_output.h5mu -t ${threshold} --cleanser
        """
}