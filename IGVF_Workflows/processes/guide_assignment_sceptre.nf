process guide_assignment_sceptre {
  
    input:
    path mudata_input

    output:
    path "guide_assignment.csv", emit: guide_assignment_matrix

    script:
    """
      assign_grnas_sceptre.R ${mudata_input}
    """
}
