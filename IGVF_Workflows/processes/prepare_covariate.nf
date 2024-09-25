
process prepare_covariate {
    
    cache 'lenient'

    input:
    val covariate_list

    output:
    path "cov_string.txt", emit: covariate_string
    path "parse_covariate.csv", emit: parsed_covariate_file
    
    script:
    def jsonString = groovy.json.JsonOutput.toJson(covariate_list).replaceAll('"', '\\\\"')
    """
    parse_covariate.py \"$jsonString\"
    prepare_formula.py parse_covariate.csv
    """

}
