#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAY_HELLO {
    output:
    stdout

    script:
    """
    echo "EPI2ME Installation erfolgreich!"
    """
}

workflow {
    SAY_HELLO()
    SAY_HELLO.out.view()
}
